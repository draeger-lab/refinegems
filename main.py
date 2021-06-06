#!/usr/bin/env python
import yaml
import refinegems as rg
import cobra

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
    
    else:
        model_cobra, errors = cobra.io.sbml.validate_sbml_model(config['model'])

        if (model_cobra != None):
            model_libsbml = rg.load_model_libsbml(config['model'])
            name, reac, metab, genes = rg.initial_analysis(model_libsbml)
            print('Name: ' + name)
            print('# reactions: ' + str(reac))
            print('# metabolites: ' + str(metab))
            print('# genes: ' + str(genes))
            #score = rg.run_memote(model_cobra)
            #print('Memote score: ' + str(score))
            if (config['media_db'] != None):
                for medium in config['media']:
                    essential, missing, growth, dt = rg.growth_simulation(model_cobra, model_libsbml, config['media_db'], medium)
                    print('---')
                    print('Growth was tested on ' + medium)
                    print('The following exchanges are needed for growth and not in the medium:')
                    print(essential)
                    print('With them the growth value is ' + str("%.2f" %growth) + ' mmol/gDW·h')
                    print('The corresponding doubling time is ' + str("%.2f" %dt)+ ' min')
                    print('The following exchanges are in the medium but not in the model:')
                    print(missing)


if __name__ == "__main__":
    main()
