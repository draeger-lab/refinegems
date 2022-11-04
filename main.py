#!/usr/bin/env python

import yaml
import refinegems as rg
import cobra
import pandas as pd
from datetime import date

__author__ = "Famke Baeuerle"

def main():
    """main function to run the program"""
    print("Report main properties of a GEM")
    print("Author:", __author__)
    today = date.today().strftime("%Y%m%d")
    
    with open('config.yaml') as f:
        config = yaml.safe_load(f)

    if (config['keggpathways']):
        non_kegg = rg.pathways.kegg_pathways(config['model'], config['kegg_path'])
        print('The following reactions have no KEGG annotation and were not added to any pathway-group: ' + str(non_kegg))
        model, errors = cobra.io.sbml.validate_sbml_model(config['kegg_path'])
        print(errors)
        
    elif (config['sboterms']):
        model_libsbml = rg.load.load_model_libsbml(config['model'])
        rg.sboann.sbo_annotation_write(model_libsbml, config['sbo_path'])
        model, errors = cobra.io.sbml.validate_sbml_model(config['sbo_path'])
        print(errors)
        
    elif (config['polish_carveme']):
        model_libsbml = rg.load.load_model_libsbml(config['model'])
        rg.polish.polish_carveme_bigg(model_libsbml, config['polish_path'], config['entrez_email'])
        model, errors = cobra.io.sbml.validate_sbml_model(config['polish_path'])
        print(errors)
        
    elif (config['charge_corr']):
        model_libsbml = rg.load.load_model_libsbml(config['model'])
        rg.charges.correct_charges_modelseed(model_libsbml, config['charge_path'], config['modelseedpath'], config['charge_report_path'])
        model, errors = cobra.io.sbml.validate_sbml_model(config['charge_path'])
        print(errors)
        
    elif(config['man_cur']):
        model_libsbml = rg.load.load_model_libsbml(config['model'])
        if config['man_cur_type'] == 'gapfill':
            gapfill = rg.load.load_manual_gapfill(config['man_cur_table'])
            model = rg.curate.add_reactions_from_table(model_libsbml, gapfill, config['entrez_email'])
            rg.load.write_to_file(model, config['man_cur_path'])
            model, errors = cobra.io.sbml.validate_sbml_model(config['man_cur_path'])
            print(errors)
        elif config['man_cur_type'] == 'metabs':
            man_ann = rg.load.load_manual_annotations(config['man_cur_table'])
            model = rg.curate.update_annotations_from_table(model_libsbml, man_ann)
            model = rg.curate.update_annotations_from_others(model)
            rg.load.write_to_file(model, config['man_cur_path'])
            model, errors = cobra.io.sbml.validate_sbml_model(config['man_cur_path'])
            print(errors)
    
    else:
        model_cobra, errors = cobra.io.sbml.validate_sbml_model(config['model'])
        print(errors)

        if (model_cobra != None):
            model_libsbml = rg.load.load_model_libsbml(config['model'])
            name, reac, metab, genes = rg.investigate.initial_analysis(model_libsbml)
            orphans, deadends, disconnected = rg.investigate.get_orphans_deadends_disconnected(model_cobra)
            mass_unbal, charge_unbal = rg.investigate.get_mass_charge_unbalanced(model_cobra)
            egc = rg.investigate.get_egc(model_cobra)
            
            if (config['memote']):
                score = rg.investigate.run_memote(model_cobra)
                
            if (config['genecomp']):
                genecomp = rg.genecomp.kegg_gene_comp(model_libsbml, config['organismid'], config['biggreactions'], config['gff_file'])
                
            if(config['modelseed']):
                charge_mismatch, formula_mismatch = rg.modelseed.compare_to_modelseed(config['modelseedpath'], model_cobra)
            
            if (config['media_db'] != None):
                growth_sim = rg.get_growth_selected_media(model_cobra, config['media_db'], config['media'])
                
            if (config['multiple']):
                growth_all = rg.comparison.simulate_all(config['multiple_paths'], config['media_db'], config['media'])
                with pd.ExcelWriter(config['out_path'] + 'growth_' + str(today) +'.xlsx') as writer:  
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
                print(egc)
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
                with pd.ExcelWriter(config['out_path'] + name + '_' + str(today) + '.xlsx') as writer:  
                    model_params.to_excel(writer, sheet_name='model params', index=False)
                    growth_sim.to_excel(writer, sheet_name='growth simulation', index=False)
                    egc.to_excel(writer, sheet_name='EGC test', index=False)
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
                model_info.to_csv(name + '_modelinfo.csv', index=False)
                growth_sim.to_csv(name +'_growthsim.csv', index=False)
                egc.to_csv(name + '_egc.csv', index=False)
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