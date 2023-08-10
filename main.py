#!/usr/bin/env python

import os
import click
import cobra
import logging
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
    
    config = rg.io.save_user_input(configpath)
    today = date.today().strftime("%Y%m%d")
    
    print('The following logs are saved to '+ config['out_path'] + 'rg_' + str(today) + '.log')

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s:%(name)s] %(message)s",
        handlers=[
            logging.FileHandler(config['out_path'] + 'rg_' + str(today) + '.log'),
            logging.StreamHandler()
        ]
    )
    logging.getLogger('cobra').setLevel(logging.WARNING)
    logging.getLogger('requests').setLevel(logging.WARNING)
    logging.getLogger('urllib3').setLevel(logging.WARNING)
    logging.info('----------- New run of refineGEMs -----------')
    
    # check if the output directory is already present, if not create it
    if not os.path.isdir(config['out_path']):
        logging.info('Given out_path is not yet a directory, creating ' + config['out_path'])
        os.makedirs(config['out_path'])
    if (config['visualize']):
        dir = os.path.join(config['out_path'] + 'visualization/')
        if not os.path.isdir(dir): 
            os.makedirs(dir)
        
    logging.info('Your output will be saved to ' + config['out_path']) 
    
    if (config['multiple']):
        logging.info('Growth simulation for multiple models: ')
        models_cobra = rg.io.load_multiple_models(config['multiple_paths'], 'cobra')
        growth_all = rg.comparison.simulate_all(models_cobra, config['media'], config['growth_basis'], config['anaerobic_growth'])
        growth_prefix = 'anaerobic_growth_' if config['anaerobic_growth'] else 'growth_'
        growth_all.to_csv(config['out_path'] + growth_prefix + str(today) + '_' + config['growth_basis'] + '.csv', index=False)
        logging.info('Multiple model growth simulation results are saved to ' +  growth_prefix + str(today) + '_' + config['growth_basis'] + '.csv')
        
        # visualizations
        if (config['visualize']):
            models_libsbml = rg.io.load_multiple_models(config['multiple_paths'], 'libsbml')
            ini_plot = rg.comparison.plot_initial_analysis(models_libsbml).get_figure()
            sbo_fig_all = rg.comparison.plot_rea_sbo_multiple(models_libsbml).get_figure()
            venn_reac = rg.comparison.plot_venn(models_cobra, 'reaction', True).get_figure()
            venn_metab = rg.comparison.plot_venn(models_cobra, 'metabolite', True).get_figure()
            heatmap = rg.comparison.plot_heatmap_dt(growth_all[['model', 'medium', 'doubling_time [min]']])
            native_heatmap = rg.comparison.plot_heatmap_native(growth_all)
            # saving them
            sbo_fig_all.savefig(config['out_path'] + 'visualization/' + 'all_ReacPerSBO_' + str(today) + '.png', bbox_inches='tight')
            venn_reac.savefig(config['out_path'] + 'visualization/' + 'all_ReacOverlap_' + str(today) + '.png', bbox_inches='tight')
            venn_metab.savefig(config['out_path'] + 'visualization/' + 'all_MetabOverlap_' + str(today) + '.png', bbox_inches='tight')
            heatmap_dt_prefix = 'heatmap_dt_additives_anaerobic_' if config['anaerobic_growth'] else 'heatmap_dt_additives_'
            heatmap.savefig(config['out_path'] + 'visualization/' + heatmap_dt_prefix + str(today) + '.png')
            native_heatmap_prefix = 'heatmap_native_anaerobic_' if config['anaerobic_growth'] else 'heatmap_native_'
            native_heatmap.savefig(config['out_path'] + 'visualization/' + native_heatmap_prefix + str(today) + '.png', bbox_inches='tight')
            ini_plot.savefig(config['out_path'] + 'visualization/' + 'model_status_' + str(today) + '.png', bbox_inches='tight')
    
    if (config['single']):        
        try:    
            model_cobra, errors = cobra.io.sbml.validate_sbml_model(config['model'])
            logging.info(errors)
        except (OSError):
            model_cobra = None
            logging.info('No or no valid model given, please enter a valid path in the model field in the config file.')

        if (config['keggpathways']):
            model_libsbml, non_kegg = rg.pathways.kegg_pathways(config['model'])
            file = open(config['out_path'] + model_libsbml.getId() + '_reac_wo_kegg_' + str(today) + '.txt','w')
            for reaction in non_kegg:
                file.write(reaction+"\n")
            file.close()
            logging.info('Kegg Pathways were added to the model as groups. Reactions that have no KEGG annotation are denoted in ' + model_libsbml.getId() + '_reac_wo_kegg.txt')

        else:
            model_libsbml = rg.io.load_model_libsbml(config['model'])
        
        if (config['sboterms']):
            if (config['visualize']):
                sbo_fig = rg.investigate.plot_rea_sbo_single(model_libsbml).get_figure()
                # saving the created visualizations
                sbo_fig.savefig(config['out_path'] + 'visualization/' + str(model_cobra.id) + '_ReacPerSBO_beforeUpdate_' + str(today) + '.png', bbox_inches='tight')
            model_libsbml = rg.sboann.sbo_annotation(model_libsbml)
            logging.info('SBO Terms updated for ' + model_libsbml.getId())
            
        if (config['charge_corr']):
            model_libsbml, multiple_charges = rg.charges.correct_charges_modelseed(model_libsbml)
            pd.DataFrame.from_dict(multiple_charges, orient='index').to_csv(config['out_path'] + model_libsbml.getId() + '_mulchar_' + str(today) + '.csv', sep=',', header=False)
            logging.info('Charges were corrected for ' + model_libsbml.getId() + '. A table with metabolites with multiple charges can be found under ' + model_libsbml.getId() + '_mulchar_' + str(today) + '.csv')
            
        if(config['man_cur']):
            if config['man_cur_type'] == 'gapfill':
                gapfill = rg.io.load_manual_gapfill(config['man_cur_table'])
                model_libsbml = rg.curate.add_reactions_from_table(model_libsbml, gapfill, config['entrez_email'])
                logging.info('Manual gap filling was done for ' + model_libsbml.getId())
            elif config['man_cur_type'] == 'metabs':
                man_ann = rg.io.load_manual_annotations(config['man_cur_table'])
                model_libsbml = rg.curate.update_annotations_from_table(model_libsbml, man_ann)
                model_libsbml = rg.curate.update_annotations_from_others(model_libsbml)
                logging.info('Manual update of annotations was done for ' + model_libsbml.getId())
                
        if config['gap_analysis'] and config['gapfill_model']:
            filename = f'{config["out_path"]}{model_libsbml.getId()}_gap_analysis_{str(today)}'
            if config['gap_analysis_params'].get('db_to_compare') not in ('BioCyc', 'KEGG+BioCyc'):
                logging.warning('Currently, only the result from the \'BioCyc\' or \'KEGG+BioCyc\' runs can be directly added to a model.')
                gap_analysis = rg.gapfill.gap_analysis(model_libsbml, config['gap_analysis_params'], filename)
                logging.info(f'Gap analysis for {model_libsbml.getId()} with {config["gap_analysis_params"].get("db_to_compare")} was performed.')
                logging.info(f'Complete Excel table is in file: {filename}.')
            else:
                gapfill = rg.gapfill.gapfill(model_libsbml, config['gap_analysis_params'], filename)
                gap_analysis_stats = gapfill[0][0]
                logging.info(f'Statistics on missing entites for {model_libsbml.getId()}:')
                logging.info(gap_analysis_stats)
                logging.info(f'Complete Excel table is in file: {filename}.')
                model_libsbml = gapfill[-1]
                logging.info(f'Gaps were filled in {model_libsbml.getId()}.')
        elif config['gap_analysis']:
            filename = f'{config["out_path"]}{model_libsbml.getId()}_gap_analysis_{str(today)}'
            gap_analysis = rg.gapfill.gap_analysis(model_libsbml, config['gap_analysis_params'], filename)
            logging.info(f'Gap analysis for {model_libsbml.getId()} with {config["gap_analysis_params"].get("db_to_compare")} was performed.')
            if  config["gap_analysis_params"].get("db_to_compare") != 'KEGG':
                logging.info(f'Statistics on missing entites for {model_libsbml.getId()}:')
                logging.info(gap_analysis[0])
            logging.info(f'Complete Excel table is in file: {filename}.')
        elif config['gapfill_model']:
            model_libsbml = rg.gapfill.gapfill_model(model_libsbml, config['gap_analysis_file'])
            logging.info(f'Gaps were filled in {model_libsbml.getId()}.')
        
        if (config['polish']):
            model_libsbml = rg.polish.polish(model_libsbml, config['entrez_email'], config['id_db'], config['gff_file'], 
                                             config['protein_fasta'], config['lab_strain'], config['organismid'], config['out_path'])
            logging.info(model_libsbml.getId() + ' has been polished')
            
        if (config['biomass']):
            result = rg.biomass.check_normalise_biomass(model_cobra)
            if result:
                model_libsbml = result
                logging.info(model_libsbml.getId() + '\'s biomass function has been checked.')
        
        mods = [config['keggpathways'], config['sboterms'], config['charge_corr'], config['man_cur'], config['gapfill_model'], config['polish'], config['biomass']]
        
        if any(mods):
            if config['model_out'] == 'stdout':   
                config['model_out'] = config['out_path'] + model_libsbml.getId() + '_modified_' + str(today) + '.xml'
                
            rg.io.write_to_file(model_libsbml, config['model_out'])
            
            if model_cobra is not None:                                          
                try:    
                    model_cobra, errors = cobra.io.sbml.validate_sbml_model(config['model_out'])
                    logging.info(errors)
                except (OSError):
                    model_cobra = None
                    logging.info('Model was invalidated during curation steps.')

        if (model_cobra != None):
            logging.info(model_cobra.id + ' will be investigated.')
            name, reac, metab, genes = rg.investigate.initial_analysis(model_libsbml)
            orphans, deadends, disconnected = rg.investigate.get_orphans_deadends_disconnected(model_cobra)
            mass_unbal, charge_unbal = rg.investigate.get_mass_charge_unbalanced(model_cobra)
            egc = rg.investigate.get_egc(model_cobra)
            if (config['visualize']):
                logging.info('All visualizations can be found in the subfolder "visualization".')
                sbo_fig = rg.investigate.plot_rea_sbo_single(model_libsbml).get_figure()
                
                # saving the created visualizations
                sbo_fig.savefig(config['out_path'] + 'visualization/' + str(model_cobra.id) + '_ReacPerSBO_' + str(today) + '.png', bbox_inches='tight')
            
            if (config['memote']):
                score = rg.investigate.get_memote_score(rg.investigate.run_memote(model_cobra))
                
            if(config['modelseed']):
                charge_mismatch, formula_mismatch = rg.modelseed.compare_to_modelseed(model_cobra)
            
            if (config['media'] != None):
                growth_sim = rg.growth.get_growth_selected_media(model_cobra, config['media'], config['growth_basis'], config['anaerobic_growth'])
            
            if (config['memote'] == True):
                information = [[name], [reac], [metab], [genes], [score], orphans, deadends, disconnected, mass_unbal, charge_unbal]
                model_params = pd.DataFrame(information, ['model name', '#reactions', '#metabolites', '#genes', 'memote score', 'orphans', 'deadends', 'disconnected', 'mass unbalanced', 'charge unbalanced']).T
            else:
                information = [[name], [reac], [metab], [genes], orphans, deadends, disconnected, mass_unbal, charge_unbal]
                model_params = pd.DataFrame(information, ['model name', '#reactions', '#metabolites', '#genes', 'orphans', 'deadends', 'disconnected', 'mass unbalanced', 'charge unbalanced']).T
            with pd.ExcelWriter(config['out_path'] + name + '_' + str(today) + '.xlsx') as writer:  
                model_params.to_excel(writer, sheet_name='model params', index=False)
                growth_sim_name = 'anaerobic growth simulation' if config['anaerobic_growth'] else 'growth simulation'
                growth_sim.to_excel(writer, sheet_name=growth_sim_name, index=False)
                egc.to_excel(writer, sheet_name='EGC test', index=False)
                if(config['modelseed']):
                    charge_mismatch.to_excel(writer, sheet_name='charge mismatches', index=False)
                    formula_mismatch.to_excel(writer, sheet_name='formula mismatches', index=False)
                logging.info('Single model growth simulation results are saved to ' + name + '_' + str(today) + '.xlsx')
        else:
            logging.info('No valid model, investigation aborted!')

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        logging.critical(e, exc_info=True)