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

    if (config['keggpathways']):
        rg.kegg_pathways(config['model'], config['kegg_path'])
        model, errors = cobra.io.sbml.validate_sbml_model(config['kegg_path'])
        print(errors)
        
    elif (config['sboterms']):
        model_libsbml = rg.load_model_libsbml(config['model'])
        rg.sbo_annotation(model_libsbml, config['database_user'], config['database_name'], config['sbo_path'])
        model, errors = cobra.io.sbml.validate_sbml_model(config['sbo_path'])
        print(errors)
        
    elif (config['polish_carveme']):
        model_libsbml = rg.load_model_libsbml(config['model'])
        rg.polish_carveme(model_libsbml, config['polish_path'], config['entrez_email'])
        model, errors = cobra.io.sbml.validate_sbml_model(config['polish_path'])
        print(errors)
        
    elif (config['charge_corr']):
        model_libsbml = rg.load_model_libsbml(config['model'])
        mulchar = rg.correct_charges(model_libsbml, config['charge_path'], config['modelseedpath'])
        model, errors = cobra.io.sbml.validate_sbml_model(config['charge_path'])
        print(errors)
        print(mulchar) # hier muss ich noch eine bessere Lösung finden, klappt aber erstmal
    
    else:
        model_cobra, errors = cobra.io.sbml.validate_sbml_model(config['model'])

        if (model_cobra != None):
            model_libsbml = rg.load_model_libsbml(config['model'])
            name, reac, metab, genes = rg.initial_analysis(model_libsbml)
            orphans, deadends, disconnected = rg.get_orphans_deadends_disconnected(model_cobra)
            mass_unbal, charge_unbal = rg.get_mass_charge_unbalanced(model_cobra)
            
            if (config['memote']):
                score = rg.run_memote(model_cobra)
                
            if (config['genecomp']):
                genecomp = rg.genecomp(model_libsbml, config['organismid'], config['biggreactions'], config['gff_file'])
                
            if(config['modelseed']):
                charge_mismatch, formula_mismatch = rg.modelseed(config['modelseedpath'], model_cobra)
            
            if (config['media_db'] != None):
                growth_sim = rg.get_growth_selected_media(model_cobra, config['media_db'], config['media'])
                
            if (config['multiple']):
                growth_all = rg.simulate_all(config['multiple_paths'], config['media_db'], config['media'])
                with pd.ExcelWriter('simulate_all.xlsx') as writer:  
                    growth_all.to_excel(writer, index=False)

            if (config['output'] == 'cl'):
                print('---')
                print('Model name: ' + name)
                print('# reactions: ' + str(reac))
                print('# metabolites: ' + str(metab))
                print('# genes: ' + str(genes))
                if (config['memote']): print('Memote score: ' + str(score))
                print('Orphan metabolites: ' + str(orphans))
                print('Deadend metabolites: ' + str(deadends))
                print('Disconnected metabolites: ' + str(disconnected))
                print('Mass unbalanced reactions: ' + str(mass_unbal))
                print('Charge unbalanced reactions: ' + str(charge_unbal))
                print(growth_sim)
                if(config['genecomp']): print(genecomp)
                if(config['modelseed']):
                    print(charge_mismatch)
                    print(formula_mismatch)
                
            if (config['output'] == 'xlsx'): # excel file
                if (config['memote'] == True):
                    information = [[name], [reac], [metab], [genes], [score], orphans, deadends, disconnected, mass_unbal, charge_unbal]
                    model_params = pd.DataFrame(information, ['model name', '#reactions', '#metabolites', '#genes', 'memote score', 'orphans', 'deadends', 'disconnected', 'mass unbalanced', 'charge unbalanced']).T
                else:
                    information = [[name], [reac], [metab], [genes], orphans, deadends, disconnected, mass_unbal, charge_unbal]
                    model_params = pd.DataFrame(information, ['model name', '#reactions', '#metabolites', '#genes', 'orphans', 'deadends', 'disconnected', 'mass unbalanced', 'charge unbalanced']).T
                with pd.ExcelWriter(name + '_refinegems.xlsx') as writer:  
                    model_params.to_excel(writer, sheet_name='model params', index=False)
                    growth_sim.to_excel(writer, sheet_name='growth simulation', index=False)
                    if(config['genecomp']):
                        genecomp.to_excel(writer, sheet_name='gene comparison', index=False)
                    if(config['modelseed']):
                        charge_mismatch.to_excel(writer, sheet_name='charge mismatches', index=False)
                        formula_mismatch.to_excel(writer, sheet_name='formula mismatches', index=False)
            
            if (config['output'] == 'csv'): # csv file
                print('---')
                print('Model name: ' + name)
                print('# reactions: ' + str(reac))
                print('# metabolites: ' + str(metab))
                print('# genes: ' + str(genes))
                if (config['memote'] == True): print('Memote score: ' + str(score))
                model_info = pd.DataFrame([orphans, deadends, disconnected, mass_unbal, charge_unbal], ['orphans', 'deadends', 'disconnected', 'mass unbalanced', 'charge unbalanced']).T
                #print('Orphan metabolites: ' + str(orphans))
                #print('Deadend metabolites: ' + str(deadends))
                #print('Disconnected metabolites: ' + str(disconnected))
                #print('Mass unbalanced reactions: ' + str(mass_unbal))
                #print('Charge unbalanced reactions: ' + str(charge_unbal))
                model_info.to_csv(name + '_modelinfo.csv', index=False)
                growth_sim.to_csv(name +'_growthsim.csv', index=False)
                if(config['genecomp']):
                    genecomp.to_csv(name +'_genecomp.csv', index=False)
                if(config['modelseed']):
                    charge_mismatch.to_csv(name + '_charge_mismatch.csv', index=False)
                    formula_mismatch.to_csv(name + '_formula_mismatch.csv', index=False)
        else:
            print(errors)
    
    print("Gem Curation Finished!")

if __name__ == "__main__":
    main()
