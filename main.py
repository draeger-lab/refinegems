#!/usr/bin/env python

import yaml
import os
import refinegems as rg
import cobra
import pandas as pd
import matplotlib.pyplot as plt
from datetime import date

__author__ = "Famke Baeuerle and Gwendolyn O. Gusak"

def main():
    """main function to run the program"""
    print("Report main properties of a GEM")
    print("Author:", __author__)
    today = date.today().strftime("%Y%m%d")
    
    rg.databases.initialise_database()
    
    with open('config.yaml') as f:
        config = yaml.safe_load(f)
    
    # check if the output directory is already present, if not create it
    dir = os.path.join(config['out_path'] + 'visualization/')
    if not os.path.isdir(config['out_path']):
        os.makedirs(config['out_path'])
    if not os.path.isdir(dir): 
        os.makedirs(dir)

    if (config['keggpathways']):
        non_kegg = rg.pathways.kegg_pathways(config['model'], config['kegg_path'])
        print('The following reactions have no KEGG annotation and were not added to any pathway-group: ' + str(non_kegg))
        model, errors = cobra.io.sbml.validate_sbml_model(config['kegg_path'])
        print(errors)
        
    elif (config['sboterms']):
        model_libsbml = rg.io.load_model_libsbml(config['model'])
        rg.sboann.sbo_annotation_write(model_libsbml, config['sbo_path'])
        model, errors = cobra.io.sbml.validate_sbml_model(config['sbo_path'])
        print(errors)
        
    elif (config['polish']):
        model_libsbml = rg.io.load_model_libsbml(config['model'])
        rg.polish.polish(model_libsbml, config['polish_path'], config['entrez_email'], config['id_db'], config['protein_fasta'], config['lab_strain'])
        model, errors = cobra.io.sbml.validate_sbml_model(config['polish_path'])
        print(errors)
        
    elif (config['charge_corr']):
        model_libsbml = rg.io.load_model_libsbml(config['model'])
        rg.charges.correct_charges_modelseed(model_libsbml, config['charge_path'], config['charge_report_path'])
        model, errors = cobra.io.sbml.validate_sbml_model(config['charge_path'])
        print(errors)
        
    elif(config['man_cur']):
        model_libsbml = rg.io.load_model_libsbml(config['model'])
        if config['man_cur_type'] == 'gapfill':
            gapfill = rg.io.load_manual_gapfill(config['man_cur_table'])
            model = rg.curate.add_reactions_from_table(model_libsbml, gapfill, config['entrez_email'])
            rg.io.write_to_file(model, config['man_cur_path'])
            model, errors = cobra.io.sbml.validate_sbml_model(config['man_cur_path'])
            print(errors)
        elif config['man_cur_type'] == 'metabs':
            man_ann = rg.io.load_manual_annotations(config['man_cur_table'])
            model = rg.curate.update_annotations_from_table(model_libsbml, man_ann)
            model = rg.curate.update_annotations_from_others(model)
            rg.io.write_to_file(model, config['man_cur_path'])
            model, errors = cobra.io.sbml.validate_sbml_model(config['man_cur_path'])
            print(errors)
    
    else:
        if (config['multiple']):
            models_cobra = rg.io.load_multiple_models(config['multiple_paths'], 'cobra')
            growth_all = rg.comparison.simulate_all(models_cobra, config['media'], config['growth_basis'])
            growth_all.to_csv(config['out_path'] + 'growth_' + str(today) + '_' + config['growth_basis'] + '.csv', index=False)
            # visualizations
            if (config['visualize']):
                models_libsbml = rg.io.load_multiple_models(config['multiple_paths'], 'libsbml')
                ini_plot = rg.comparison.plot_initial_analysis(models_libsbml).get_figure()
                sbo_fig_all = rg.comparison.plot_rea_sbo_multiple(models_libsbml).get_figure()
                venn_reac = rg.comparison.plot_venn(models_cobra, 'reaction', True).get_figure()
                venn_metab = rg.comparison.plot_venn(models_cobra, 'metabolite', True).get_figure()
                heatmap = rg.comparison.plot_heatmap_dt(growth_all[['model', 'medium', 'doubling_time [min]']])
                binary_heatmap = rg.comparison.plot_heatmap_binary(growth_all)
                # saving them
                sbo_fig_all.savefig(config['out_path'] + 'visualization/' + 'all_ReacPerSBO_' + str(today) + '.png', bbox_inches='tight')
                venn_reac.savefig(config['out_path'] + 'visualization/' + 'all_ReacOverlap_' + str(today) + '.png', bbox_inches='tight')
                venn_metab.savefig(config['out_path'] + 'visualization/' + 'all_MetabOverlap_' + str(today) + '.png', bbox_inches='tight')
                heatmap.savefig(config['out_path'] + 'visualization/' + 'heatmap_dt_additives_' + str(today) + '.png')
                binary_heatmap.savefig(config['out_path'] + 'visualization/' + 'heatmap_native_' + str(today) + '.png', bbox_inches='tight')
                ini_plot.savefig(config['out_path'] + 'visualization/' + 'model_status_' + str(today) + '.png', bbox_inches='tight')
                                                    
        try:    
            model_cobra, errors = cobra.io.sbml.validate_sbml_model(config['model'])
            print(errors)
        except (OSError):
            model_cobra = None
            print('Either no or no valid model given, please enter a valid path in the model field in the config file.')

        if (model_cobra != None):
            model_libsbml = rg.io.load_model_libsbml(config['model'])
            name, reac, metab, genes = rg.investigate.initial_analysis(model_libsbml)
            orphans, deadends, disconnected = rg.investigate.get_orphans_deadends_disconnected(model_cobra)
            mass_unbal, charge_unbal = rg.investigate.get_mass_charge_unbalanced(model_cobra)
            egc = rg.investigate.get_egc(model_cobra)
            sbo_fig = rg.investigate.plot_rea_sbo_single(model_libsbml).get_figure()
            
            # saving the created visualizations
            sbo_fig.savefig(config['out_path'] + 'visualization/' + str(model_cobra.id) + '_ReacPerSBO_' + str(today) + '.png', bbox_inches='tight')
            
            if (config['memote']):
                score = rg.investigate.run_memote(model_cobra)
                
            if (config['gapfill_analysis'] and config['gapfill_model']):
                filename = config['out_path'] + name + '_gapfill_analysis_' + str(today) + '.xlsx'
                gapfill_complete = rg.gapfill.gapfill(model_libsbml, config['gapfill_analysis_params'], filename) 
                gapfill_analysis = gapfill_complete[:-1]
                model = gapfill_complete[-1]
                rg.io.write_to_file(model, config['gapfill_model_out'])
                model, errors = cobra.io.sbml.validate_sbml_model(config['gapfill_model_out'])
                print(errors)   
            elif (config['gapfill_analysis']):
                filename = config['out_path'] + name + '_gapfill_analysis_' + str(today) + '.xlsx'
                gapfill_analysis = rg.gapfill.gapfill_analysis(model_libsbml, config['gapfill_analysis_params'], filename)
            elif (config['gapfill_model']):
                model = rg.gapfill.gapfill_model(model_libsbml, config['gapfill_model_in'])
                rg.io.write_to_file(model, config['gapfill_model_out'])
                model, errors = cobra.io.sbml.validate_sbml_model(config['gapfill_model_out'])
                print(errors)   
                
            if(config['modelseed']):
                charge_mismatch, formula_mismatch = rg.modelseed.compare_to_modelseed(model_cobra)
            
            if (config['media'] != None):
                growth_sim = rg.growth.get_growth_selected_media(model_cobra, config['media'], config['growth_basis'])

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
                if(config['gapfill_analysis']) or (config['gapfill_analysis'] and config['gapfill_model']): 
                    if type(gapfill_analysis) == tuple:
                        print('BioCyc - Statistics on missing entities:')
                        print(gapfill_analysis[0])
                        if len(gapfill_analysis) == 6:
                            print(f'Complete Excel table is in file: {config["out_path"]}{name}_gapfill_analysis_{str(today)}.xlsx')
                        else:
                            print(f'Complete Excel table is in file: {config["out_path"]}{name}_gapfill_analysis_{str(today)}.xlsx')
                    else:
                        print(f'Complete Excel table is in file: {config["out_path"]}{name}_gapfill_analysis_{str(today)}.xlsx')
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
                if(config['modelseed']):
                    charge_mismatch.to_csv(name + '_charge_mismatch.csv', index=False)
                    formula_mismatch.to_csv(name + '_formula_mismatch.csv', index=False)
    
    print("Gem Curation Finished!")

if __name__ == "__main__":
    main()