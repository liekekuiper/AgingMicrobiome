import os
import pandas as pd
import numpy as np
import statsmodels.api as sm
from lifelines import CoxPHFitter
from lifelines.utils import k_fold_cross_validation
from microbiome_utils_fixed_alex3 import *
from patsy import dmatrices

# Fill in the right filepaths
label = '16s/metagenomics'
taxonomy = '2022.10.taxonomy.asv.tsv.qza'
tree = '2022.10.phylogeny.asv.nwk.qza'
feature_table = 'feature.table.gg2-2022.10.qza'
metadata = 'metadata.txt'
threads = int(6)
modsfile = 'mods_agingmicrobiome.txt'
outsfile = 'outs_agingmicrobiome.txt'
cohortname = 'cohort.txt'
subsfile = 'subsets_agingmicrobiome.txt'
factors_str = 'NA'

#----------------------------------------------------#
##   Do not change anything underneath these lines  ##
#----------------------------------------------------#

if factors_str.startswith("'") and factors_str.endswith("'"):
    factors_str = factors_str[1:-1]

factors = factors_str.split(',')

print(threads)
print(factors)
print('imported environments')

results_df = pd.DataFrame(columns=[
    'Datasplit', 'Outcome', 'Variable', 'Model', 'N', 'Ncases', 'Coefficient',
    'Std.Error', 'HR', 'LL','UL','t.value', 'P'
])

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

for subs in subsets:
    print(f"Performing analyses of {subs}")
    for out in outcomes:
        for model in models:
            original_model = model
            if out not in ["age", "mortality"]:
                model = f"{model}+age"
            if "men" in subs:
                model = model.replace("sex+", "").replace("sex", "")
            model_terms = model.split('+')
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
                label=label
            )
            datafile.columns = datafile.columns.str.replace('-', '_')
            datafile.columns = datafile.columns.str.replace(' ', '_')
            datafile.columns = datafile.columns.str.replace('[^0-9a-zA-Z_]', '', regex=True)

            print('created dataset, now continue with analyses')
            variables = [
                col for col in datafile.columns 
                if col not in variables_store and not any(col.startswith(term + '_') for term in model_terms)
            ]
            for var in variables:
                if out != "mortality":
                    formula = f"{out} ~ {model} + {var}"
                    datafile[out] = pd.to_numeric(datafile[out])
                    y, X = dmatrices(formula, data=datafile, return_type='dataframe')
                    lmmodel = sm.OLS(y, X).fit()
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
                    results_df = pd.concat([results_df, result_row], ignore_index=True)
                else:
                    datafile['followup'] = datafile['studytime'] + datafile['age']
                    covariates = model_terms + [var]
                    covariates = [cov.strip() for cov in covariates]
                    covariates_use = [col for col in datafile.columns if any(cov in col for cov in covariates)]
                    cph = CoxPHFitter()
                    try:
                        cph.fit(datafile[covariates_use + ['followup', 'mortality']], duration_col='followup', event_col='mortality')
                        coefficients = cph.summary
                        variable = var
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
                        print(f"Failed CoxPH fit for variable '{var}' in subset '{subs}' with outcome '{out}'. Error: {e}")
                        print(f"Covariates used: {covariates_use}")
                        continue
            results_df.to_csv(f"./intermediatefiles/Results_{label}_{subs}_{cohort}_{pd.Timestamp.today().date()}.csv", index=False)
            model = original_model
    print(f"Finished analyses of {subs}")

results_df.to_csv(f"Results_{cohort}_{label}_{pd.Timestamp.today().date()}.csv", index=False)

if __name__ == '__main__':
    process()
