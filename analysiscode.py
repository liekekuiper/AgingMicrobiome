import pandas as pd
import numpy as np
import statsmodels.api as sm
from lifelines import CoxPHFitter
from lifelines.utils import k_fold_cross_validation
from microbiome_utils import *

# Prepare dataframe
results_df = pd.DataFrame(columns=[
    'Datasplit', 'Outcome', 'Variable', 'Model', 'N', 'Ncases', 'Coefficient',
    'Std.Error', 'HR', 't.value', 'P'
])

print(f"Performing analyses of {subs}")
for out in outcomes:
    for var in variables:
        for model in models:
            if out not in ["age", "mortality"]:
                model = f"{model}+age"
            if "men" in subs:
                model = model.replace("sex+", "").replace("sex", "")
            datafile = process(taxonomy, tree, feature_table, output, threads, metadata, model, sub)
            if out != "mortality":
                formula = f"{out} ~ {model} + {var}"
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
                    't.value': np.nan,
                    'P': coefficients.loc[variable, 'p']
                }])

            results_df = pd.concat([results_df, result_row], ignore_index=True)

print(f"Finished analyses of {subs}")
results_df.to_csv(f"Results_{subs}_{cohort}_{pd.Timestamp.today().date()}.csv", index=False)

