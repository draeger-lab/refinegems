#!/usr/bin/env python
from refinegems.genecomp import genecomp
import yaml
import refinegems as rg
import cobra
import pandas as pd

__author__ = "Famke Baeuerle"

def main():
    """main function to run the program"""
    print("Report main properties of a GEM")
    print("Author:", __author__)
    
    with open('config.yaml') as f:
        config = yaml.safe_load(f)

    if (config['keggpathways'] != None):
        rg.kegg_pathways(config['model'], config['keggpathways'])
        model, errors = cobra.io.sbml.validate_sbml_model(config['keggpathways'])
        print(errors)
        
    elif (config['sboterms']):
        model_libsbml = rg.load_model_libsbml(config['model'])
        rg.sbo_annotation(model_libsbml, config['database_user'], config['database_name'], config['new_filename'])
        model, errors = cobra.io.sbml.validate_sbml_model(config['new_filename'])
        print(errors)
    
    else:
        model_cobra, errors = cobra.io.sbml.validate_sbml_model(config['model'])

        if (model_cobra != None):
            model_libsbml = rg.load_model_libsbml(config['model'])
            name, reac, metab, genes = rg.initial_analysis(model_libsbml)
            
            if (config['memote']):
                score = rg.run_memote(model_cobra)
                
            if (config['genecomp']):
                genecomp = rg.genecomp(model_libsbml, config['organismid'], config['biggreactions'], config['gff_file'])
                
            if (config['media_db'] != None):
                df_list = []
                for medium in config['media']: # ACHTUNG, so funktioniert das natürlich nicht :D
                    model_cobra = rg.load_model_cobra(config['model'])
                    model_libsbml = rg.load_model_libsbml(config['model'])
                    essential, missing, growth, dt = rg.growth_simulation(model_cobra, model_libsbml, 'media/media_db.csv', medium)
                    exchanges = [[medium], essential, missing, [growth], [dt]]
                    df_temp = pd.DataFrame(exchanges, ['name', 'essential', 'missing', 'growth_value [mmol/gDW·h]', 'doubling_time [min]']).T
                    df_list.append(df_temp)
                growth_sim = pd.concat(df_list)

            if (config['output'] == 'cl'):
                print('---')
                print('Model name: ' + name)
                print('# reactions: ' + str(reac))
                print('# metabolites: ' + str(metab))
                print('# genes: ' + str(genes))
                if (config['memote']): print('Memote score: ' + str(score))
                print(growth_sim)
                if(config['genecomp']): print(genecomp)
                
            if (config['output'] == 'xlsx'): # excel file
                if (config['memote'] == True):
                    information = [[name], [reac], [metab], [genes], [score]]
                    model_params = pd.DataFrame(information, ['model name', '#reactions', '#metabolites', '#genes', 'memote score']).T
                else:
                    information = [[name], [reac], [metab], [genes]]
                    model_params = pd.DataFrame(information, ['model name', '#reactions', '#metabolites', '#genes']).T
                with pd.ExcelWriter(name + '_refinegems.xlsx') as writer:  
                    model_params.to_excel(writer, sheet_name='model params', index=False)
                    growth_sim.to_excel(writer, sheet_name='growth simulation', index=False)
                    if(config['genecomp']):
                        genecomp.to_excel(writer, sheet_name='gene comparison', index=False)
            
            if (config['output'] == 'csv'): # csv file
                print('---')
                print('Model name: ' + name)
                print('# reactions: ' + str(reac))
                print('# metabolites: ' + str(metab))
                print('# genes: ' + str(genes))
                if (config['memote'] == True):
                    print('Memote score: ' + str(score))
                growth_sim.to_csv(name +'_growthsim.csv', index=False)
                if(config['genecomp']):
                    genecomp.to_csv(name +'_genecomp.csv', index=False)
    
    print("Gem Curation Finished!")

if __name__ == "__main__":
    main()
