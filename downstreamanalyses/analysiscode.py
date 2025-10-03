import os
import pandas as pd
import numpy as np
import statsmodels.api as sm
from lifelines import CoxPHFitter
from lifelines.utils import k_fold_cross_validation  # kept even if not used, safe to remove if you prefer
from microbiome_utils import *  # uses process() and helpers
from patsy import dmatrices

# reads paths and settings from environment variables provided by sbatch --export
taxonomy = os.getenv('taxonomy')  # path to taxonomy artifact
tree = os.getenv('tree')  # path to tree artifact
feature_table = os.getenv('feature_table')  # path to feature table artifact
metadata = os.getenv('metadata')  # path to metadata tsv
threads = int(os.getenv('threads', '6'))  # number of threads, defaults to 6 if not set
modsfile = os.getenv('modsfile')  # file with model formulas (one per line)
outsfile = os.getenv('outsfile')  # file with outcomes (one per line)
cohortname = os.getenv('cohortname')  # file with cohort label (single line)
subsfile = os.getenv('subsfile')  # file with subset names (one per line)
label = os.getenv('label', '16s')  # gets label from environment, defaults to '16s'
factors_str = os.getenv('factors', 'NA')  # study-specific categorical covariates as comma-separated string

# removes outer single quotes if they were passed in sbatch as "'A,B'"
if factors_str.startswith("'") and factors_str.endswith("'"):
    factors_str = factors_str[1:-1]

factors = factors_str.split(',')  # turns "A,B" into ["A","B"]

print(threads)  # quick log of threads
print(factors)  # quick log of factors
print('imported environments')

# prepares an empty results dataframe with consistent columns
results_df = pd.DataFrame(columns=[
    'Datasplit', 'Outcome', 'Variable', 'Model', 'N', 'Ncases', 'Coefficient',
    'Std.Error', 'HR', 'LL', 'UL', 't.value', 'P'
])

# loads model list, outcome list, cohort label, and subset list
with open(modsfile, 'r') as file:
    models = file.read().splitlines()  # each line is a base model string like 'sex+ppump'
with open(outsfile, 'r') as file:
    outcomes = file.read().splitlines()  # each line is an outcome name
with open(cohortname, 'r') as file:
    cohort = file.readline().strip()  # single cohort label
with open(subsfile, 'r') as file:
    subsets = file.read().splitlines()  # each line is a subset like 'all' or 'men'

print('read files')

# stores original metadata columns to help exclude them from the features to test
variables_store = pd.read_csv(metadata, sep='\t', dtype=str).columns.str.lower()

print('stored variables, will start for loop')

# loops over all subset × outcome × model combinations
for subs in subsets:
    print(f"Performing analyses of {subs}")
    for out in outcomes:
        for model in models:
            original_model = model  # keeps the original for resetting later

            # adds age unless the outcome itself is age or mortality
            if out not in ["age", "mortality"]:
                model = f"{model}+age"

            # removes sex from the model when analyzing only men
            if "men" in subs:
                model = model.replace("sex+", "").replace("sex", "")

            model_terms = model.split('+')  # splits model into individual terms

            # builds the analysis dataset using QIIME2 artifacts and diversity metrics
            datafile = process(
                taxonomy=taxonomy,
                tree=tree,
                feature_table=feature_table,
                output=cohort,
                threads=threads,
                metadata=metadata,
                model=model,
                sub=subs,
                out=out,
                factors=factors,
                label=label  # passes label so species-level is included for metagenomics
            )

            # standardizes column names to avoid patsy/formula issues
            datafile.columns = datafile.columns.str.replace('-', '_')  # replaces dashes with underscores
            datafile.columns = datafile.columns.str.replace(' ', '_')  # replaces spaces with underscores
            datafile.columns = datafile.columns.str.replace('[^0-9a-zA-Z_]', '', regex=True)  # removes other special chars

            print('created dataset, now continue with analyses')

            # selects candidate variables to test (excludes original metadata cols and model covariate dummies)
            variables = [
                col for col in datafile.columns
                if col not in variables_store and not any(col.startswith(term + '_') for term in model_terms)
            ]

            # runs OLS for continuous outcomes and Cox PH for mortality
            for var in variables:
                if out != "mortality":
                    # builds formula "out ~ model + var"
                    formula = f"{out} ~ {model} + {var}"
                    datafile[out] = pd.to_numeric(datafile[out])  # ensures numeric outcome
                    y, X = dmatrices(formula, data=datafile, return_type='dataframe')  # constructs design matrices
                    lmmodel = sm.OLS(y, X).fit()  # fits linear model
                    coefficients = lmmodel.summary2().tables[1]  # extracts coefficient table
                    variable = var  # the tested feature

                    # appends one row with estimates for the tested feature
                    result_row = pd.DataFrame([{
                        'Datasplit': subs,
                        'Outcome': out,
                        'Variable': var,
                        'Model': model,
                        'N': len(lmmodel.fittedvalues),
                        'Ncases': np.nan,
                        'Coefficient': coefficients.loc[variable, 'Coef.'],
                        'Std.Error': coefficients.loc[variable, 'Std.Err.'],
                        'HR': np.nan,
                        'LL': coefficients.loc[variable, '[0.025'],
                        'UL': coefficients.loc[variable, '0.975]'],
                        't.value': coefficients.loc[variable, 't'],
                        'P': coefficients.loc[variable, 'P>|t|']
                    }])
                    results_df = pd.concat([results_df, result_row], ignore_index=True)
                else:
                    # creates attained-age time scale: followup = studytime + age
                    datafile['followup'] = datafile['studytime'] + datafile['age']

                    # builds covariate list including model terms and the tested feature
                    covariates = model_terms + [var]
                    covariates = [cov.strip() for cov in covariates]  # trims whitespace

                    # maps to the actual encoded columns (incl. one-hot dummies)
                    covariates_use = [col for col in datafile.columns if any(cov in col for cov in covariates)]

                    cph = CoxPHFitter()
                    try:
                        # fits Cox model using attained-age follow-up and mortality as event
                        cph.fit(
                            datafile[covariates_use + ['followup', 'mortality']],
                            duration_col='followup',
                            event_col='mortality'
                        )
                        coefficients = cph.summary  # coefficient table from lifelines
                        variable = var  # the tested feature

                        # appends one row with estimates for the tested feature
                        result_row = pd.DataFrame([{
                            'Datasplit': subs,
                            'Outcome': out,
                            'Variable': var,
                            'Model': model,
                            'N': cph._n_examples,
                            'Ncases': cph.event_observed.sum(),
                            'Coefficient': coefficients.loc[variable, 'coef'],
                            'Std.Error': coefficients.loc[variable, 'se(coef)'],
                            'HR': coefficients.loc[variable, 'exp(coef)'],
                            'LL': coefficients.loc[variable, 'exp(coef) lower 95%'],
                            'UL': coefficients.loc[variable, 'exp(coef) upper 95%'],
                            't.value': np.nan,
                            'P': coefficients.loc[variable, 'p']
                        }])
                        results_df = pd.concat([results_df, result_row], ignore_index=True)
                    except Exception as e:
                        print(f"Failed CoxPH fit for variable '{var}' in subset '{subs}' with outcome '{out}'. Error: {e}")  # logs failure reason
                        print(f"Covariates used: {covariates_use}")  # logs which columns were used
                        continue  # moves on to the next variable

            # writes an intermediate CSV after each model to help monitor progress
            results_df.to_csv(
                f"./intermediatefiles/Results_{label}_{subs}_{cohort}_{pd.Timestamp.today().date()}.csv",
                index=False
            )

            model = original_model  # resets model string for the next loop iteration

    print(f"Finished analyses of {subs}")

# writes the final combined results CSV at the end
results_df.to_csv(f"Results_{cohort}_{label}_{pd.Timestamp.today().date()}.csv", index=False)
