# This file is used for developement of the refinegems module
# %%
import pandas as pd
import refinegems as rg
import cobra

model_cobra, errors = cobra.io.sbml.validate_sbml_model('models/CStr_20210518.xml')

if (model_cobra != None):
    model_libsbml = rg.load_model_libsbml('models/CStr_20210518.xml')
    name, reac, metab, genes = rg.initial_analysis(model_libsbml)
    print('Name: ' + name)
    print('# reactions: ' + str(reac))
    print('# metabolites: ' + str(metab))
    print('# genes: ' + str(genes))
    
information = [[name], [reac], [metab], [genes]]
df1 = pd.DataFrame(information, ['model name', '#reactions', '#metabolites', '#genes']).T
#df1.to_excel('test.xlsx', sheet_name='model parameters')

# %%
df_list = []
for medium in ['SNM3', 'M9']:
    model_cobra, errors = cobra.io.sbml.validate_sbml_model('models/CStr_20210518.xml')
    model_libsbml = rg.load_model_libsbml('models/CStr_20210518.xml')
    essential, missing, growth, dt = rg.growth_simulation(model_cobra, model_libsbml, 'media/media_db.csv', medium)
    exchanges = [[medium], essential, missing, [growth], [dt]]
    df_temp = pd.DataFrame(exchanges, ['name', 'essential', 'missing', 'growth_value', 'doubling_time']).T
    df_list.append(df_temp)

df2 = pd.concat(df_list)
#df.to_excel('test.xlsx', sheet_name='growth simulation')

#%%
print(df1)
df2

df2.to_csv('test.csv', index=False)

# %%
pd.read_csv('test.csv')

# %%
with pd.ExcelWriter('test.xlsx') as writer:  
    df1.to_excel(writer, sheet_name='model params', index=False)
    df2.to_excel(writer, sheet_name='growth simulation', index=False)
