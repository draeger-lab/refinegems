"""Entry points to the code from the command line.
"""

__author__ = 'Carolin Brune, Gwendolyn O. Döbel'

################################################################################
# Requirements
################################################################################

from traitlets import default
from typing import Union
import refinegems as rg
import click
from pathlib import Path

################################################################################
# Entry points
################################################################################

@click.group()
@click.help_option("--help", "-h")
@click.version_option()
def cli():
   """refineGEMs - A toolbox for fast curation and analysis of constraint-based metabolic models
   
   This tool provides scripts and functions to curate and analyse a strain-specific GEM with information from multiple 
   databases. In addition, it offers functions for analysis of multiple GEMs at the same time as well as comparison 
   reports.
   """

# ------
# Set-up
# ------
@cli.group()
def setup():
   """Set-up tools, folder structure and more for running the program.
   """

# Get a config file
# -----------------
@setup.command()
@click.option('--filename', '-f', default='config.yaml', type=str, 
              show_default=True, help='Name (path) to save the config file under.')
@click.option('--type', '-t', default='media', type=click.Choice(['media', 'refinegems']), 
              help='Type of config file to download. Either a file for media configuration or a file to run a refinement pipeline.')
def config(filename,type):
    """Download a configuration file (.yaml).

    Download a configuration file to edit for running the complete
    workflow or provide the media configuration.
    """
    rg.utility.set_up.download_config(filename, type)

# @TODO: Download a file with all media for a specific database?
# Handle medium / media database
# ------------------------------
@setup.command()
@click.option('--list', is_flag=True, default=False, help='List the names of the media in the database.')
#@click.option('--copy', is_flag=False, flag_value='./media.csv', default=None, type=click.Path(exists=False), 
# help='Produce and save a copy of the media database.')
def medium(list): #,copy
    """Access the in-build media database.

    Can be used to either check the available media or to make a copy for further use.
    """
    if list:
        possible_media = rg.utility.io.load_a_table_from_database('medium', False)['name'].to_list()
        media = '|'.join(possible_media)
        print(media)
'''
    if copy:
        db = specimen.classes.medium.load_media_db()
        specimen.classes.medium.save_db(db,click.format_filename(copy))
'''

# build a pan-core model
# ----------------------
# @TEST
@setup.command()
@click.argument('models',nargs=-1,type=click.Path())
@click.option('-o','--based-on', required=False, type=click.Choice(['id']), default='id',help='Option on how to combine the models.')
@click.option('-n','--name', required=False, type=str, default='pan-core-model',help='Name of the new pan-core model.')
@click.option('-g','--keep-genes',is_flag=True, default=False)
@click.option('--rcomp', '--resolve-compartments',is_flag=True)
@click.option('-d', '--dir', required=False, type=click.Path(), default='', help='Path to the output dir.')
def build_pancore(models, based_on, name, keep_genes, rcomp,dir):
   """Build a pan-core model.
   """
   pancore_mod = rg.analysis.core_pan.generate_core_pan_model(models, based_on, name, not keep_genes, rcomp)
   rg.utility.io.write_model_to_file(pancore_mod,Path(dir,name +'.xml'))

# --------------
# Polish a model
# --------------
@cli.group()
def polish():
	"""Polish a model created by an automatic reconstruction pipeline. CLeans up the notes fields, changes qualifiers and 
 	annotations to be MIRIAM-compliant and fixes the initial biomass equation to add up to 1mmol/gDW/h.
	"""

# Polish an automatically generated reconstruction
# ------------------------------------------------
@polish.command()
@click.argument('model', type=str)
@click.argument('email', type=str)
@click.argument('path', type=str)
@click.option('-i','--id_db', default='BiGG', type=str, help='Main database where identifiers in model come from')
@click.option('-r', '--refseq_gff', default=None, type=str, help='Path to RefSeq GFF file of organism')
@click.option('-p', '--protein_fasta', default=None, type=str, help='File used as input for CarveMe')
@click.option('-l', '--lab_strain', default=False, type=bool, help='True if the strain was sequenced in a local lab')
@click.option('-k', '--kegg_organism_id', default=None, type=str, help='KEGG organism identifier')
def polish(model,email,path,id_db,refseq_gff,protein_fasta,lab_strain,kegg_organism_id):
   """Completes all steps to polish a model

   (Tested for models having either BiGG or VMH identifiers.)
   """
   rg.curation.polish.polish(model, email, id_db, refseq_gff, protein_fasta, lab_strain, kegg_organism_id, path)


# ----------------------------------------------------------------------------------
# Find and fill gaps in a model automatically/Fill gaps with manually created tables
# ----------------------------------------------------------------------------------
@cli.group()
def gaps():
   """If requested finds gaps in a model based on the gene products in the model and either files from biocyc or via the
   KEGG API.
   If requested fills gaps either via the automatically obtained missing entities or via a manually created file 
   provided by the user.
   """

# Find gaps via genes
# -------------------
@gaps.command()
@click.argument('modelpath', type=str)
@click.argument('gff_file', type=str)
@click.argument('organismid', type=str)
@click.argument('gapfill_params', type=dict)
@click.argument('filename', type=str)
def find(modelpath,gff_file,organismid,gapfill_params,filename):
   """Find gaps in a model based on the genes/gene products of the underlying organism
   """
   model = rg.utility.io.load_model(modelpath, 'libsbml')
   rg.curation.gapfill.gap_analysis(model, gff_file, organismid, gapfill_params, filename)

# Fill gaps via file
# ------------------
@gaps.command()
@click.argument('model', type=str)
@click.argument('gap_analysis_results', type=Union[str, tuple])
def fill(modelpath,gap_analysis_result):
   """Fill gaps in a model based on a user-provided input file
   """
   model = rg.utility.io.load_model(modelpath, 'libsbml')
   rg.curation.gapfill.gapfill_model(model, gap_analysis_result)

# Find and fill gaps via genes
# ----------------------------
@gaps.command()
@click.argument('modelpath', type=str)
@click.argument('gapfill_params', type=dict)
@click.argument('filename', type=str)
def autofill(modelpath,gapfill_params,filename):
   """Automatically find and fill gaps based on the genes/gene products
   """
   model = rg.utility.io.load_model(modelpath, 'libsbml')
   rg.curation.gapfill.gapfill(model, gapfill_params, filename)


# --------------
# Refine a model
# --------------
@cli.group()
def refine():
   """Refine a model. Includes steps like biomass, charges, SBO annotation, reaction direction correction and addition 
   of Pathways and further gene product annotations.
   """
   
@refine.command()
@click.argument('modelpath', type=str)
@click.option('-c', '--cycles', default=10, type=int, help='Maximal number of optiomisation cycles that will be run.')
def biomass(modelpath,cycles):
   """Changes the biomass reaction to be consistent
   """
   model = rg.utility.io.load_model(modelpath, 'cobra')
   rg.curation.biomass.check_normalise_biomass(model, cycles)
   
@refine.command()
@click.argument('modelpath', type=str)
def charges(modelpath):
   """Changes the charges in a model by comparison to the ModelSEED database
   """
   model = rg.utility.io.load_model(modelpath, 'libsbml')
   rg.curation.charges.correct_charges_modelseed(model)

# @TODO
@refine.command()
@click.argument('modelpath', type=str)
def direction(modelpath,data):
   """Checks & if necessary, corrects the direction for each reaction in a model
   """
   model = rg.utility.io.load_model(modelpath, 'cobra')
   rg.curation.polish.check_direction(model, data)

# @IDEA maybe more fitting in annot group?
@refine.command()
@click.argument('modelpath', type=str)
def pathways(modelpath):
   """Add KEGG pathways to a model
   """
   rg.curation.pathways.kegg_pathways(modelpath)

# Annotation-related clean-up & Additional annotations
# ---------------------------
@refine.group()
def annot():
	"""Clean-up annotations by adding MIRIAM-compliant qualifiers, changing URIs into the correct form, adding more 
 	specific SBO annotations and synchronising all annotations to BioCyc. Additionally, adds gene product annotations 
  	either from the GFF file of the organism or via KEGG.
	"""
 
@annot.command()
@click.argument('modelpath', type=str)
def sboterms(modelpath):
   """Calls SBOannotator to enhance the SBO terms in a model
   """
   model = rg.utility.io.load_model(modelpath, 'libsbml')
   rg.sboann.sbo_annotation(model)


# -----------------------------------------------
# Analyse a model/multiple models including a comparison
# -----------------------------------------------
@cli.group()
def analyse():
	"""Analyse a model by testing for auxotrophies and growth on different media along with finding EGCs and looking at 
 	overall model statistics. 
	"""

@analyse.command()
@click.argument('modelpath', type=str)
@click.option('-s', '--score-only', is_flag=True, default=False, help='Specifies if memote is only run to return the score')
@click.option('-f','--file', required=False, type=str, default='memote.html', help='Name or path to save the output to. Only relevent if -s it not set.')
def memote(modelpath,score_only,file):
   """Perform a memote analysis.
   """
   model = rg.utility.io.load_model(modelpath, 'cobra')
   if score_only:
      rg.analysis.investigate.get_memote_score(rg.analysis.investigate.run_memote(model,type='json',return_res=True))
   else:
      rg.analysis.investigate.run_memote(model,save_res=file)
 
 # @TODO add colour option
@analyse.command()
@click.argument('modelpath', type=str)
@click.option('-d', '--dir', required=False, type=click.Path(), default='', help='Path to the output dir.')
def pathways(modelpath,dir):
   """Analysis of pathways contained in a model
   """
   model = rg.utility.io.load_model(modelpath, 'cobra')
   report = rg.curation.pathways.kegg_pathway_analysis(model)
   report.save(dir)

# @TODO: accept multiple models
@analyse.command()
@click.argument('modelpath', type=click.Path(exists=True))
@click.option('-d', '--dir', required=False, type=click.Path(), default='', help='Path to the output dir.')
@click.option('-c', '--colors', required=False, type=str, default='YlGn', help='Abbreviation of a matplotlib colour palette.')
def stats(modelpath,dir,colors):
   """Generate a report on the statistics of a model.
   """
   model = rg.utility.io.load_model(modelpath,'cobra')
   report = rg.classes.reports.ModelInfoReport(model)
   report.save(dir,colors)

# @TEST
@analyse.command()
@click.argument('modelpath', type=click.Path(exists=True))
@click.argument('pcpath', type=click.Path(exists=True))
@click.option('-b','--based-on', type=click.Choice(['id']),required=False, default='id',help='Option on how to compare the models.')
@click.option('-d', '--dir', required=False, type=click.Path(), default='', help='Path to the output dir.')
def pancore(modelpath, pcpath, based_on,dir):
   """Compare a model to a pan-core model.
   """
   model = rg.utility.io.load_model(modelpath,'cobra')
   pcmodel = rg.utility.io.load_model(pcpath,'cobra')
   report = rg.analysis.core_pan.compare_to_core_pan(model,pcmodel,based_on)
   report.save(dir)

# analyse growth
# --------------
@analyse.group()
def growth():
   """Analyse the growth under different conditions."""

# @TEST
# -> does it need to run without a media config as well?
@growth.command()
@click.argument('modelpaths', nargs=-1, type=click.Path(exists=True))
@click.option('-m','--media',required=True,type=click.Path(exists=True), help='Path to a media config file.')
@click.option('-n', '--namespace', required=False, type=click.Choice(['BiGG']), default='BiGG', help='Namespace to use for the model.')
@click.option('-d', '--dir', required=False, type=click.Path(), default='', help='Path to the output dir.')
@click.option('-c', '--colors', required=False, type=str, default='YlGn', help='Abbreviation of a matplotlib colour palette.')
def simulate(modelpaths,media,namespace,dir,colors):
   """Simulate the growth of the given model vs. media.
   """
   report = rg.analysis.growth.growth_analysis(modelpaths, media, namespace,retrieve='report')
   report.save(to=dir, color_palette=colors)

# @TEST
@growth.command()
@click.argument('modelpath', type=click.Path(exists=True))
@click.option('-m','--media',required=True,type=click.Path(exists=True), help='Path to a media config file.')
@click.option('-n', '--namespace', required=False, type=click.Choice(['BiGG']), default='BiGG', help='Namespace to use for the model.')
@click.option('-d', '--dir', required=False, type=click.Path(), default='', help='Path to the output dir.')
@click.option('-c', '--colors', required=False, type=str, default='YlGn', help='Abbreviation of a matplotlib colour palette.')
def auxotrophies(modelpath,media,namespace,dir,colors):
   """Test for auxotrophies for the 20 proteinogenic amino acids.
   """
   medialist, supps = rg.analysis.growth.read_media_config(media)
   model = rg.utility.io.load_model(modelpath, 'cobra')
   report = rg.analysis.growth.test_auxotrophies(model,medialist,namespace)
   report.save(dir,colors)

@growth.command()
@click.argument('modelpath', type=click.Path(exists=True))
@click.option('-e','--element', type=str, required=True, help='Element to perform the source test with. Needs to be a valid chemical, elemental Symbol.')
@click.option('-s','--substances', type=str, required=False, default=None, multiple=True, help='Add substances from the database to test against. If none are given, resrs against the while database.')
@click.option('-m', '--medium', required=False, type=str, default=None, help='Name of a medium from the database to use for the testing. If not given, uses the one from the model.')
@click.option('-n', '--namespace', required=False, type=click.Choice(['BiGG']), default='BiGG', help='Namespace to use for the model.')
@click.option('-d', '--dir', required=False, type=click.Path(), default='', help='Path to the output dir.')
@click.option('-c', '--colors', required=False, type=str, default='YlGn', help='Abbreviation of a matplotlib colour palette.')
def sources(modelpath, element, substances, medium, namespace, dir, colors):
   """Test growth on different sources for an element.
   """ 
   model = rg.utility.io.load_model(modelpath,'cobra')
   report = rg.analysis.growth.test_growth_with_source(model,element,substances,medium,namespace)
   report.save(dir,color_palette=colors)


@growth.command()
@click.argument('modelpath', type=click.Path(exists=True))
@click.option('-o','--objective', required=False, type=click.Choice(['flux','medium','exchanges']), default='flux',help='Set the type of minimal medium to be calculated: minimal fluxes, minimal compounds or based on the model exchange reactions.')
@click.option('-r', '--growth-rate', required=False, type=float, default=0.5, help='Minimal rate to be reached by the minimal medium. The smaller the rate the more expensive the calculation.')
@click.option('-d', '--dir', required=False, type=click.Path(), default='', help='Path to the output dir.')
def minimal_medium(modelpath,objective, growth_rate, dir):
   """Calculate the minimal medium of a model.

   Either set the fluxes minimal, the number of compounds based on current medium or the number of exchange reactions.
   """
    
   model = rg.utility.io.load_model(modelpath,'cobra')
   medium = rg.analysis.growth.model_minimal_medium(model, objective, growth_rate)
   medium.export_to_file(dir=dir)