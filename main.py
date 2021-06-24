#!/usr/bin/env python
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
        
    elif (config['sboterms'] != None):
        model_libsbml = rg.load_model_libsbml(config['model'])
        rg.sbo_annotation(model_libsbml, config['sboterms'][0], config['sboterms'][1], config['sboterms'][2])
        model, errors = cobra.io.sbml.validate_sbml_model(config['sboterms'][2])
        print(errors)
    
    else:
        model_cobra, errors = cobra.io.sbml.validate_sbml_model(config['model'])

        if (model_cobra != None):
            model_libsbml = rg.load_model_libsbml(config['model'])
            name, reac, metab, genes = rg.initial_analysis(model_libsbml)
            
            if (config['memote'] == True):
                score = rg.run_memote(model_cobra)
                
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

            if (config['output'][0] == 'command_line'):
                print('---')
                print('Model name: ' + name)
                print('# reactions: ' + str(reac))
                print('# metabolites: ' + str(metab))
                print('# genes: ' + str(genes))
                if (config['memote'] == True):
                    print('Memote score: ' + str(score))
                print(growth_sim)
                
            if (config['output'][0] != 'command_line'):
                if (config['output'][1] == 1): # excel file
                    if (config['memote'] == True):
                        information = [[name], [reac], [metab], [genes], [score]]
                        model_params = pd.DataFrame(information, ['model name', '#reactions', '#metabolites', '#genes', 'memote score']).T
                    else:
                        information = [[name], [reac], [metab], [genes]]
                        model_params = pd.DataFrame(information, ['model name', '#reactions', '#metabolites', '#genes']).T
                    with pd.ExcelWriter(config['output'][0]) as writer:  
                        model_params.to_excel(writer, sheet_name='model params', index=False)
                        growth_sim.to_excel(writer, sheet_name='growth simulation', index=False)
                if (config['output'][1] == 2): # csv file
                    print('---')
                    print('Model name: ' + name)
                    print('# reactions: ' + str(reac))
                    print('# metabolites: ' + str(metab))
                    print('# genes: ' + str(genes))
                    if (config['memote'] == True):
                        print('Memote score: ' + str(score))
                    growth_sim.to_csv(config['output'][0], index=False)
    
    print("Gem Curation Finished!")

if __name__ == "__main__":
    main()
