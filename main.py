#!/usr/bin/env python

import os
import click
import cobra
import refinegems as rg
import pandas as pd
import matplotlib.pyplot as plt
from datetime import date

__author__ = "Famke Baeuerle and Gwendolyn O. Gusak"

@click.command()
@click.option('-c', '--configpath', required=True, prompt='Enter path to config file or press Enter if you want to create one.',
              help='Path to file containing configurations to run refineGEMs. An example file can be found in config.yaml.', default='ENTER')

def main(configpath=None):
    """main function to run the program"""
    print("RefineGEMs provides functions to curate and investigate genome-scale metabolic models!")
    print("Author:", __author__)
    today = date.today().strftime("%Y%m%d")
    
   # click.echo('Config file loaded from ' + click.format_filename(configpath))
    config = rg.io.save_user_input(configpath)
    rg.databases.initialise_database()
    
    # check if the output directory is already present, if not create it
    dir = os.path.join(config['out_path'] + 'visualization/')
    if not os.path.isdir(config['out_path']):
        print('Given out_path is not yet a directory, creating ' + config['out_path'])
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
            if (config['visualize']):
                sbo_fig = rg.investigate.plot_rea_sbo_single(model_libsbml).get_figure()
                
                # saving the created visualizations
                sbo_fig.savefig(config['out_path'] + 'visualization/' + str(model_cobra.id) + '_ReacPerSBO_' + str(today) + '.png', bbox_inches='tight')
            
            if (config['memote']):
                score = rg.investigate.run_memote(model_cobra)
                
            if (config['gapfill_analysis'] and config['gapfill_model']):
                gapfill_analysis = rg.gapfill.gapfill(model_libsbml, config['gapfill_analysis_params'])    
            elif (config['gapfill_analysis']):
                gapfill_analysis = rg.gapfill.gapfill_analysis(model_libsbml, config['gapfill_analysis_params'])
            elif (config['gapfill_model']):
                rg.gapfill.gapfill_model(model_libsbml, config['gapfill_model_file'])
                
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
                        print('BioCyc - Missing genes for reactions table:')
                        print(gapfill_analysis[1])
                        print('BioCyc - Missing metabolites with BiGG ID table:')
                        print(gapfill_analysis[2])
                        print('BioCyc - Missing metabolites without BiGG ID table:')
                        print(gapfill_analysis[3])
                        print('BioCyc - Missing reactions with BiGG ID table:')
                        print(gapfill_analysis[4])
                        if len(gapfill_analysis) == 6:
                            print('KEGG - Missing reactions with BiGG ID table:')
                            print(gapfill_analysis[6])
                    else:
                        print('KEGG - Gap fill analysis')
                        print(gapfill_analysis)
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
                    if(config['gapfill_analysis']) and type(gapfill_analysis) != tuple:
                        gapfill_analysis.to_excel(writer, sheet_name='KEGG gap fill analysis', index=False)
                    if(config['modelseed']):
                        charge_mismatch.to_excel(writer, sheet_name='charge mismatches', index=False)
                        formula_mismatch.to_excel(writer, sheet_name='formula mismatches', index=False)
                if(config['gapfill_analysis']) or (config['gapfill_analysis'] and config['gapfill_model']): 
                    if type(gapfill_analysis) == tuple:
                        with pd.ExcelWriter(config['out_path'] + name + '_gapfill_analysis_' + str(today) + '.xlsx') as writer:
                            gapfill_analysis[0].to_excel(writer, sheet_name='gap fill statistics', index=False)
                            gapfill_analysis[1].to_excel(writer, sheet_name='genes', index=False)
                            gapfill_analysis[2].to_excel(writer, sheet_name='metabolites', index=False)
                            gapfill_analysis[3].to_excel(writer, sheet_name='metabolites without BiGG IDs', index=False)
                            gapfill_analysis[4].to_excel(writer, sheet_name='reactions', index=False)
                            if len(gapfill_analysis) == 6:
                                gapfill_analysis[5].to_excel(writer, sheet_name='KEGG reactions', index=False)
            
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
                if(config['gapfill_analysis']) or (config['gapfill_analysis'] and config['gapfill_model']):
                    if type(gapfill_analysis) == tuple:
                        gapfill_analysis[0].to_csv(name +'_BioCyc_analysis_statistics.csv', index=False)
                        gapfill_analysis[1].to_csv(name +'_BioCyc_analysis_genes.csv', index=False)
                        gapfill_analysis[2].to_csv(name +'_BioCyc_analysis_metabolites.csv', index=False)
                        gapfill_analysis[3].to_csv(name +'_BioCyc_analysis_metabolites_wo_BiGG.csv', index=False)
                        gapfill_analysis[4].to_csv(name +'_BioCyc_analysis_reactions.csv', index=False) 
                        if len(gapfill_analysis) == 6:
                            gapfill_analysis[5].to_csv(name +'_KEGG_gapfill_analysis.csv', index=False) 
                    else:
                        gapfill.to_csv(name +'_KEGG_gapfill_analysis.csv', index=False)
                if(config['modelseed']):
                    charge_mismatch.to_csv(name + '_charge_mismatch.csv', index=False)
                    formula_mismatch.to_csv(name + '_formula_mismatch.csv', index=False)
    
    print("Gem Curation Finished!")

if __name__ == "__main__":
    main()