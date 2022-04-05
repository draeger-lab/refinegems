# This file is used for developement of the refinegems module
# %%
from numpy import NaN
import pandas as pd
import refinegems as rg
import cobra
import re
import numpy as np
from libsbml import *
#%%

def return_modelseed_compounds(path):
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

df_modelseed = return_modelseed_compounds('modelseed/modelseed_compounds.tsv')#.to_csv('modelseed/modelseed_compounds_condensed.csv', index=False)

#%%

model = cobra.io.read_sbml_model('models/cstr_ma.xml')

#%%
charges = {}
for metab in model.metabolites:
    charges[metab.id[:-2]] = [metab.charge, metab.formula]
    
df_charges = pd.DataFrame.from_dict(charges, orient='index', columns=['charge_model', 'formula_model']).reset_index().rename(columns={'index':'BiGG'})

df_charges

#%%
df_modelseed = df_modelseed.loc[:,['charge', 'BiGG', 'formula']].rename(columns={'charge':'charge_modelseed', 'formula':'formula_modelseed'})

df_modelseed
#%%
lst_col='BiGG'
x = df_modelseed.assign(**{lst_col:df_modelseed[lst_col].str.split(';')})
df_ms = pd.DataFrame({col:np.repeat(x[col].values, x[lst_col].str.len()) for col in x.columns.difference([lst_col])}).assign(**{lst_col:np.concatenate(x[lst_col].values)})[x.columns.tolist()]
df_ms
#%%
df_comp = pd.merge(df_charges, df_ms, on='BiGG', how='left')
df_comp

#%%
def f(x):
    return True if float(x['charge_model']) == x['charge_modelseed'] else False

def g(x):
    return True if x['formula_model'] == x['formula_modelseed'] else False

df_comp['charge_match'] = df_comp.apply(f, axis=1)
df_comp['formula_match'] = df_comp.apply(g, axis=1)

#df_comp.loc[~df_comp['charge_match']].dropna(subset=['charge_modelseed'])
df_comp

#%% new stuff
model_libsbml = rg.load_model_libsbml('../Nextcloud/master_thesis/models/Cstr_17_sbo.xml')

spe = model_libsbml.getListOfSpecies()

com = pd.read_csv('modelseed/modelseed_compounds_condensed.csv', sep=',')
com

# get charges for previously uncharged metabolites
mulchar = dict()
for i in spe:
    if not i.getPlugin('fbc').isSetCharge(): # we are only interested in metab without charge
        bigg = i.getId()[2:-2]
        if len(com[com['BiGG'] == bigg]['charge'].array) == 1: #eindeutig
            charge = com[com['BiGG'] == bigg]['charge'].array[0]
            i.getPlugin('fbc').setCharge(int(charge))
        elif len(com[com['BiGG'] == bigg]['charge'].array) > 1:
            charges = com[com['BiGG'] == bigg]['charge'].array
            if all(x==charges[0] for x in charges):
                charge = charges[0]
                i.getPlugin('fbc').setCharge(int(charge))
            else:
                mulchar[bigg] = charges

print(mulchar)
            

#%%     
new_document = model_libsbml.getSBMLDocument()
writeSBMLToFile(new_document, 'Cstr.xml')
