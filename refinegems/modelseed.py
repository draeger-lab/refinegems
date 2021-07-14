import pandas as pd
import cobra
import re
import numpy as np

def get_modelseed_compounds(path):
    com = pd.read_csv(path, sep='\t')

    def get_bigg_ids(aliases):
        try: 
            aliases_list = aliases.split('|')
            bigg = [x[6:] for x in aliases_list if re.search('BiGG: .*', x)]
            return bigg[0]
        except (IndexError, AttributeError):
            return None

    com['BiGG'] = com.apply(lambda row: get_bigg_ids(row['aliases']), axis=1)

    return com.loc[:,['id', 'name', 'formula', 'mass', 'charge', 'BiGG']].dropna(subset=['BiGG'])

def get_model_charges(model):
    charges = {}
    for metab in model.metabolites:
        charges[metab.id[:-2]] = [metab.charge, metab.formula]
        
    df_charges = pd.DataFrame.from_dict(charges, orient='index', columns=['charge_model', 'formula_model']).reset_index().rename(columns={'index':'BiGG'})

    return df_charges

def get_modelseed_charges(modelseed_compounds):
    modelseed_compounds = modelseed_compounds.loc[:,['charge', 'BiGG', 'formula']].rename(columns={'charge':'charge_modelseed', 'formula':'formula_modelseed'})
    lst_col='BiGG'
    x = modelseed_compounds.assign(**{lst_col:modelseed_compounds[lst_col].str.split(';')})
    df_ms = pd.DataFrame({col:np.repeat(x[col].values, x[lst_col].str.len()) for col in x.columns.difference([lst_col])}).assign(**{lst_col:np.concatenate(x[lst_col].values)})[x.columns.tolist()]
    return df_ms

def compare_model_modelseed(model_charges, modelseed_charges):
    df_comp = pd.merge(model_charges, modelseed_charges, on='BiGG', how='left')

    def f(x):
        return True if float(x['charge_model']) == x['charge_modelseed'] else False

    def g(x):
        return True if x['formula_model'] == x['formula_modelseed'] else False

    df_comp['charge_match'] = df_comp.apply(f, axis=1)
    df_comp['formula_match'] = df_comp.apply(g, axis=1)

    return df_comp

def get_charge_mismatch(df_comp):
    return df_comp.loc[~df_comp['charge_match']].dropna(subset=['charge_modelseed'])

def get_formula_mismatch(df_comp):
    return df_comp.loc[~df_comp['formula_match']].dropna(subset=['formula_modelseed'])

def modelseed(path, model):
    ms_comp = get_modelseed_compounds(path)
    model_charges = get_model_charges(model)
    modelseed_charges = get_modelseed_charges(ms_comp)
    df_comp = compare_model_modelseed(model_charges, modelseed_charges)
    charge_mismatch = get_charge_mismatch(df_comp)
    formula_mismatch = get_formula_mismatch(df_comp)
    return charge_mismatch, formula_mismatch