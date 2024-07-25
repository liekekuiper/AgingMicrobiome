import os
import pandas as pd
import numpy as np
import statsmodels.api as sm
from lifelines import CoxPHFitter
from lifelines.utils import k_fold_cross_validation
from microbiome_utils import *
from patsy import dmatrices

# Retrieve environment variables
taxonomy = os.getenv('taxonomy')
tree = os.getenv('tree')
feature_table = os.getenv('feature_table')
metadata = os.getenv('metadata')
threadsget =  os.getenv('threads')
threads = int(threadsget)
modsfile = os.getenv('modsfile')
outsfile = os.getenv('outsfile')
cohortname = os.getenv('cohortname')
subsfile = os.getenv('subsfile')

print(threads)

print('imported environments')

# Prepare dataframe
results_df = pd.DataFrame(columns=[
    'Datasplit', 'Outcome', 'Variable', 'Model', 'N', 'Ncases', 'Coefficient',
    'Std.Error', 'HR', 'LL','UL','t.value', 'P'
])

# Load data
with open(modsfile, 'r') as file:
    models = file.read().splitlines()
with open(outsfile, 'r') as file:
    outcomes = file.read().splitlines()
with open(cohortname, 'r') as file:
    cohort = file.readline().strip()
with open(subsfile, 'r') as file:
    subsets = file.read().splitlines()

print('read files')

variables_store = pd.read_csv(metadata, sep='\t', dtype=str).columns.str.lower()

print('stored variables, will start for loop')

#Run analyses
for subs in subsets:
    print(f"Performing analyses of {subs}")
    for out in outcomes:
        for model in models:
            original_model = model
            if out not in ["age", "mortality"]:
                model = f"{model}+age"
            if "men" in subs:
                model = model.replace("sex+", "").replace("sex", "")
            datafile = process(
                taxonomy=taxonomy,
                tree=tree,
                feature_table=feature_table,
                output=cohort,
                threads=threads,
                metadata=metadata,
                model=model,
                sub=subs
            )
            print('created diversity matrices, now continue with analyses')
            variables = [col for col in datafile.columns if col not in variables_store]
            for var in variables:
                if out != "mortality":
                    formula = f"{out} ~ {model} + {var}"
                    datafile[out] = pd.to_numeric(datafile[out])
                    y, X = dmatrices(formula, data=datafile, return_type='dataframe')
                    lmmodel = sm.OLS(y, X).fit()
                    print(formula)                   
                    coefficients = lmmodel.summary2().tables[1]
                    variable = var
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
                else:
                    datafile['followup'] = datafile['studytime'] + datafile['age']
                    formula = f"followup ~ {model} + {var}"
                    cph = CoxPHFitter()
                    cph.fit(datafile, duration_col='followup', event_col='mortality')
                    coefficients = cph.summary
                    variable = var
                    result_row = pd.DataFrame([{
                        'Datasplit': subs,
                        'Outcome': out,
                        'Variable': var,
                        'Model': model,
                        'N': cph._n_examples,
                        'Ncases': cph._n_events,
                        'Coefficient': coefficients.loc[variable, 'coef'],
                        'Std.Error': coefficients.loc[variable, 'se(coef)'],
                        'HR': coefficients.loc[variable, 'exp(coef)'],
                        'LL': coefficients.loc[variable, 'exp(coef) lower 95%'],
                        'UL': coefficients.loc[variable, 'exp(coef) upper 95%'],
                        't.value': np.nan,
                        'P': coefficients.loc[variable, 'p']
                    }])
                results_df = pd.concat([results_df, result_row], ignore_index=True)
            model = original_model  # Reset the model to its original state for the next iteration
    print(f"Finished analyses of {subs}")

results_df.to_csv(f"Results_{subs}_{cohort}_{pd.Timestamp.today().date()}.csv", index=False)

if __name__ == '__main__':
    process()  # Ensuring the click context for process() to work if this script is run standalone
