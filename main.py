#!/usr/bin/env python
import click
import refinegems as rg
import cobra

__author__ = "Famke Baeuerle"

@click.command()
@click.option('-i', '--input', default='models/e_coli_core.xml', help='Path to GEM to be investigated')
@click.option('-o', '--output', default=False, help='Determine if output file should be created, default FALSE')
@click.option('-k', '--keggpathways', default=None, help='new filename/path for model with KEGG pathways')
@click.option('-m', '--medium', default='media/media_db.csv', help='Simulate growth on medium defined in csv file, default SNM3')
def main(input, output, keggpathways, medium):
    """main function to run the program"""
    print("Report main properties of a GEM")
    print("Author:", __author__)

    if (keggpathways != None):
        rg.kegg_pathways(input, keggpathways)
        model, errors = cobra.io.sbml.validate_sbml_model(keggpathways)
    
    else:
        model_cobra, errors = cobra.io.sbml.validate_sbml_model(input)

        if (model_cobra != None):
            model_libsbml = rg.load_model_libsbml(input)
            name, reac, metab, genes = rg.initial_analysis(model_libsbml)
            print('Name: ' + name)
            print('# reactions: ' + str(reac))
            print('# metabolites: ' + str(metab))
            print('# genes: ' + str(genes))
            score = rg.run_memote(model_cobra)
            print('Memote score: ' + str(score))
            if (medium != None):
                mediumname = 'SNM3'
                essential, missing, growth, dt = rg.growth_simulation(model_cobra, model_libsbml, medium, mediumname)
                print('Growth was tested on ' + mediumname)
                print('The following exchanges are needed for growth and not in the medium:')
                print(essential)
                print('With them the growth value is ' + str("%.2f" %growth) + ' mmol/gDW·h')
                print('The corresponding doubling time is ' + str("%.2f" %dt)+ ' min')
                print('The following exchanges are in the medium but not in the model:')
                print(missing)


if __name__ == "__main__":
    main()
