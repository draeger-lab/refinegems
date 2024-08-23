"""Entry points to the code from the command line.
"""

__author__ = 'Carolin Brune, Gwendolyn O. Döbel'

################################################################################
# Requirements
################################################################################

from pathlib import Path
from typing import Union, Literal

import click
import cloup
import pandas as pd
import warnings

import refinegems as rg
from refinegems.utility.io import write_model_to_file

################################################################################
# Entry points
################################################################################

@cloup.group()
@cloup.help_option("--help", "-h")
@cloup.version_option()
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

# download important stuff
# ------------------------
@setup.command()
@click.option('--filename', '-f', default='config.yaml', type=str, 
              show_default=True, help='Name (path) to save the config file under.')
@click.option('--type', '-t', show_default=True, default='media', type=click.Choice(['media', 'refinegems']), 
              help='Type of config file to download. Either a file for media configuration or a file to run a refinement pipeline.')
def config(filename,type):
    """Download a configuration file (.yaml).

    Download a configuration file to edit for running a specific function
    or provide the media configuration file.
    """
    rg.utility.set_up.download_config(filename, type)
    

@setup.command()
@click.argument('downloadtype', type=click.Choice(['SwissProt_gapfill']))
@click.option('-d', '--dir', required=False, type=click.Path(), show_default=True,  default='', help='Path to the output dir.')
@click.option('-c', '--chunksize', required=False, type=int, show_default=True, default=10, help='Chunksizre for the download in kB. Defaults to 10.')
def data(downloadtype, dir, chunksize):
   """Download file(s) needed for a given functionality of the toolbox.
   """
   dltype = downloadtype.replace('_',' ')
   rg.utility.set_up.download_url(dltype, dir,chunksize)
   

# build a pan-core model
# ----------------------
# @TEST
@setup.command()
@click.argument('models',nargs=-1,type=click.Path())
@click.option('-o','--based-on', required=False, type=click.Choice(['id']), show_default=True, default='id',help='Option on how to combine the models.')
@click.option('-n','--name', required=False, type=str, show_default=True, default='pan-core-model',help='Name of the new pan-core model.')
@click.option('-g','--keep-genes',is_flag=True, show_default=True, default=False, help='If set, the genes are kept in the pan-core model, otherwise they are deleted.')
@click.option('--rcomp', '--resolve-compartments',is_flag=True, help='If set, tries to standardise the compartment names to the c,p,e,... namespace.')
@click.option('-d', '--dir', required=False, type=click.Path(), show_default=True,  default='', help='Path to the output dir.')
def build_pancore(models, based_on, name, keep_genes, rcomp,dir):
   """Build a pan-core model.
   """
   pancore_mod = rg.analysis.core_pan.generate_core_pan_model(models, based_on, name, not keep_genes, rcomp)
   rg.utility.io.write_model_to_file(pancore_mod,Path(dir,name +'.xml'))

# ------------
# get examples
# ------------
# @DISCUSSION explain how to do manual gapfilling


# --------------------
# work on the database
# --------------------

@cli.group()
def database():
   """Access, curate, etc. the in-build database.
   """
   
@database.command()
def initialise():
   """Initialise or update the database.
   """
   rg.utility.databases.initialise_database()
   
@database.command()
@click.argument('databasename', type=click.Choice(['MetaNetX']))
@click.option('-c','--chunksize',show_default=True, default=1, help='Chunksize for the download in kB. Defaults to 1.')
def add_namespace(databasename, chunksize):
   """Add or update table for a given namespace/database to the database.
   """
   match databasename:
      case 'MetaNetX':
         rg.utility.databases.update_mnx_namespaces(chunksize=chunksize)
      case _:
         mes = f'Unknown database name: {databasename}'
         raise ValueError(mes)
         
   
@database.command()
def reset():
   """Reset the database by removing additionally downloaded tables.
   """
   rg.utility.databases.reset_database()

# -------------------
# all about the media 
# -------------------
# @TODO more functionalities / entry points

@cli.group()
def media():
   """Access the media part of the database.
   """

# @TODO: Download a file with all media for a specific database?
# Handle medium / media database
# ------------------------------
@media.command()
@click.option('--list', is_flag=True, show_default=True, default=False, help='List the names of the media in the database.')
#@click.option('--copy', is_flag=False, flag_value='./media.csv', default=None, type=click.Path(exists=False), 
# help='Produce and save a copy of the media database.')
def info(list): #,copy
    """Access information about the in-build media database.

    Can be used to either check 
    - the available media 
    -  @TODO to make a copy for further use.
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
@click.option('-i','--id-db', show_default=True, default='BiGG', type=str, help='Main database where identifiers in model come from')
@click.option('-r', '--refseq-gff', show_default=True, default=None, type=str, help='Path to RefSeq GFF file of organism') # @DISCUSSION only RefSeq or does GenBank work as well?
@click.option('-p', '--protein-fasta', show_default=True, default=None, type=str, help='File used as input for CarveMe')
@click.option('-l', '--lab-strain', show_default=True, default=False, type=bool, help='True if the strain was sequenced in a local lab')
@click.option('-k', '--kegg-organism-id', show_default=True, default=None, type=str, help='KEGG organism identifier')
def run(model,email,path,id_db,refseq_gff,protein_fasta,lab_strain,kegg_organism_id):
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
def get_gap_analysis_input(db_to_compare: Literal['KEGG', 'BioCyc']) -> dict:
   """Form to retrieve input files for :py:func:`~refinegems.curation.gapfill.gap_analysis`

   Args:
       - db_to_compare (Literal['KEGG', 'BioCyc']): 
           Database to compare model content to

   Returns:
       dict: 
           Input dictionary for :py:func:`~refinegems.curation.gapfill.gap_analysis`
   """

   parameters2inputs = {'organismid': None, 'gff_file': None, 'biocyc_files': None}

   if db_to_compare == 'KEGG' or db_to_compare == 'KEGG+BioCyc':
       parameters2inputs['organismid'] = click.prompt('Enter the KEGG Organism ID', type=str)
       parameters2inputs['gff_file'] = click.prompt('Enter the path to your organisms RefSeq GFF file', type=click.Path(exists=True))
   if db_to_compare == 'BioCyc' or db_to_compare == 'KEGG+BioCyc':
       Path0 = click.prompt('Enter the path to your BioCyc TXT file containing a SmartTable with the columns \'Accession-2\' and \'Reaction of gene\'', type=click.Path(exists=True))
       Path1 = click.prompt('Enter the path to your BioCyc TXT file containing a SmartTable with all reaction relevant information', type=click.Path(exists=True))
       Path2 = click.prompt('Enter the path to your Biocyc TXT file containing a SmartTable with all metabolite relevant information', type=click.Path(exists=True))
       Path3 = click.prompt('Enter path to protein FASTA file used as input for CarveMe', type=click.Path(exists=True))
       parameters2inputs['biocyc_files'] = [Path0, Path1, Path2, Path3]
   
   return parameters2inputs

@gaps.command()
@cloup.argument('modelpath', type=click.Path(exists=True), help='Path to model file')
@cloup.argument('filename', type=str, help='Path to output file for gap analysis result')
@cloup.option_group(
    'Reference database options',
    cloup.option('-k', '--kegg', type=bool, help='Specifies if the KEGG API should be used for the gap analysis'),
    cloup.option('-b', '--biocyc', 
                 type=bool, help='Specifies if user-provided BioCyc files should be used for the gap analysis'),
    constraint=cloup.constraints.RequireAtLeast(1)
)
def find(modelpath,kegg,biocyc,filename):
   """Find gaps in a model based on the genes/gene products of the underlying organism
   """
   
   db_to_compare = ''
   if kegg and biocyc: db_to_compare = 'KEGG+BioCyc'
   elif kegg: db_to_compare = 'KEGG'
   elif biocyc: db_to_compare = 'BioCyc'
   
   params2ins = get_gap_analysis_input(db_to_compare)
   
   model = rg.utility.io.load_model(modelpath, 'libsbml')
   rg.curation.gapfill.gap_analysis(model, db_to_compare, params2ins['gff_file'], params2ins['organismid'], params2ins['biocyc_files'], filename)

# Fill gaps via file
# ------------------
@gaps.command()
@click.argument('modelpath', type=click.Path(exists=True))
@click.argument('gap_analysis_results', type=Union[str, tuple])
def fill(modelpath,gap_analysis_result):
   """Fill gaps in a model based on a user-provided input file
   """
   model = rg.utility.io.load_model(modelpath, 'libsbml')
   rg.curation.gapfill.gapfill_model(model, gap_analysis_result)

# Find and fill gaps via genes
# ----------------------------
@gaps.command()
@cloup.argument('modelpath', type=click.Path(exists=True), help='Path to model file')
@cloup.argument('filename', type=str, help='Path to output file for gap analysis result')
@cloup.option_group(
    'Reference database options',
    cloup.option('-k', '--kegg', type=bool, help='Specifies if the KEGG API should be used for the gap analysis'),
    cloup.option('-b', '--biocyc', 
                 type=bool, help='Specifies if user-provided BioCyc files should be used for the gap analysis'),
    constraint=cloup.constraints.RequireAtLeast(1)
)
def autofill(modelpath,kegg,biocyc,filename):
   """Automatically find and fill gaps based on the genes/gene products
   """
   
   db_to_compare = ''
   if kegg and biocyc: db_to_compare = 'KEGG+BioCyc'
   elif kegg: db_to_compare = 'KEGG'
   elif biocyc: db_to_compare = 'BioCyc'
   
   params2ins = get_gap_analysis_input(db_to_compare)
   
   model = rg.utility.io.load_model(modelpath, 'libsbml')
   rg.curation.gapfill.gapfill(model, db_to_compare, params2ins['gff_file'], params2ins['organismid'], params2ins['biocyc_files'], filename)


# --------------
# Refine a model
# --------------
@cli.group()
def refine():
   """Refine a model. Includes steps like biomass, charges, SBO annotation, reaction direction correction and addition 
   of Pathways and further gene product annotations.
   """

@refine.command()
@click.argument('modelpath', type=click.Path(exists=True))
@click.option('-c', '--cycles', show_default=True, default=10, type=int, help='Maximal number of optimisation cycles that will be run.')
@click.option('--outfile','-o',required=False, type=click.Path(), show_default=True, default='biomass_corrected.xml', help='Path to save the corrected model under.')
def biomass(modelpath,cycles,outfile):
   """Changes the biomass reaction to be consistent
   """
   model = rg.utility.io.load_model(modelpath, 'cobra')
   corrected = rg.curation.biomass.check_normalise_biomass(model, cycles)
   if corrected:
      rg.utility.io.write_model_to_file(model,outfile)
   

# @TEST
@refine.command()
@click.argument('modelpath', type=click.Path(exists=True))
@click.option('--dir','-d',required=False, type=click.Path(), show_default=True, default='', help='Directory to save the output to.')
def charges(modelpath,dir):
   """Changes the charges in a model by comparison to the ModelSEED database
   """
   model = rg.utility.io.load_model(modelpath, 'libsbml')
   corr,charg = rg.curation.charges.correct_charges_modelseed(model)
   rg.utility.io.write_model_to_file(corr,Path(dir,'corrected_model.xml'))
   charg_tab = pd.DataFrame(charg)
   charg_tab.to_csv(Path(dir,'multiple_charges.csv',sep=';'))
       

# @TODO function unfinished
@refine.command()
@click.argument('modelpath', type=click.Path(exists=True))
def direction(modelpath,data):
   """Checks & if necessary, corrects the direction for each reaction in a model
   """
   model = rg.utility.io.load_model(modelpath, 'cobra')
   rg.curation.polish.check_direction(model, data)
   warnings.warn('Function not yet fully implemented - please contact developers.')
   pass


# egcs
#@TODO: Default for compartments -> Why c,p and not c,e,p?
@refine.command()
@click.argument('modelpath', type=click.Path(exists=True))
@click.option('--solver', '-s',required=False,type=click.Choice(['greedy']), show_default=True, default=None, multiple=False, help='Type of solver for the EGCs.')
@click.option('--namespace','-n',required=False, type=click.Choice(['BiGG']), show_default=True, default='BiGG',multiple=False, help='Namespace of the model.')
@click.option('--compartment','-c', required=False, type=str, show_default=True, default='c,p', help='Compartments to check, separated by comma only.')
@click.option('--outfile','-o', required=False, type=click.Path(), show_default=True, default='fixed_egcs.xml',help='Path to save edited model to, if solver has been set.')
def egcs(modelpath, solver, namespace, compartment, outfile):
   """Identify and optionally solve EGCs.
   """
   model = rg.utility.io.load_model(modelpath, 'cobra')
   match solver:
      case 'greedy':
            solver = rg.classes.egcs.GreedyEGCSolver()
            solution = solver.solve_egcs(model,namespace,compartment)
            for k,v in solution.items():
               print(k + ': ' + v)
            rg.utility.io.write_model_to_file(model,outfile)
      case _:
           solver = rg.classes.egcs.EGCSolver()
           egcs,vals = solver.find_egcs(model,True,namespace,compartment)
           print(f'EGCs: \n{str(egcs)}')
           print(f'Objective values: \n{str(vals)}')



# Annotation-related clean-up 
# & Additional annotations
# ---------------------------
@refine.group()
def annot():
	"""Add annotations to your model.
	"""
 
@annot.command()
@click.argument('modelpath', type=click.Path(exists=True))
def sboterms(modelpath):
   """Calls SBOannotator to enhance the SBO terms in a model
   """
   model = rg.utility.io.load_model(modelpath, 'libsbml')
   rg.utility.connections.run_SBOannotator(model)


@annot.command()
@click.argument('modelpath', type=click.Path(exists=True))
@click.option('-d','--dir', required=False, default='./', help='Path to the output dorectory')
def pathways(modelpath):
   """Add KEGG pathways as groups to a model
   """
   model, missing = rg.curation.pathways.kegg_pathways(modelpath) 
   with open(Path(dir, 'reac_wo_kegg_pathway_groups.txt'), 'w') as outfile:
      # save reactions with missing groups
      for line in missing:
            outfile.write(f'{line}\n')
   # save model 
   write_model_to_file(Path(dir,'model_with_added_KeggPathwayGroups.xml'))


# -----------------------------------------------
# Analyse a model/multiple models including a comparison
# -----------------------------------------------
@cli.group()
def analyse():
	"""Analyse a model by testing for auxotrophies and growth on different media along with finding EGCs and looking at 
 	overall model statistics. 
	"""

@analyse.command()
@click.argument('modelpath', type=click.Path(exists=True))
@click.option('-s', '--score-only', is_flag=True, show_default=True, default=False, help='Specifies if memote is only run to return the score')
@click.option('-f','--file', required=False, type=str, show_default=True, default='memote.html', help='Name or path to save the output to. Only relevant if -s is not set.')
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
@click.argument('modelpath', type=click.Path(exists=True))
@click.option('-d', '--dir', required=False, type=click.Path(), show_default=True, default='', help='Path to the output dir.')
def pathways(modelpath,dir):
   """Analysis of pathways contained in a model
   """
   model = rg.utility.io.load_model(modelpath, 'cobra')
   report = rg.curation.pathways.kegg_pathway_analysis(model)
   report.save(dir)


# @TODO: accept multiple models
@analyse.command()
@click.argument('modelpath', type=click.Path(exists=True))
@click.option('-d', '--dir', required=False, type=click.Path(), show_default=True, default='', help='Path to the output dir.')
@click.option('-c', '--colors', required=False, type=str, show_default=True, default='YlGn', help='Abbreviation of a matplotlib colour palette.')
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
@click.option('-b','--based-on', type=click.Choice(['id']),required=False, show_default=True, default='id',help='Option on how to compare the models.')
@click.option('-d', '--dir', required=False, type=click.Path(), show_default=True, default='', help='Path to the output dir.')
def pancore(modelpath, pcpath, based_on,dir):
   """Compare a model to a pan-core model.
   """
   model = rg.utility.io.load_model(modelpath,'cobra')
   pcmodel = rg.utility.io.load_model(pcpath,'cobra')
   report = rg.analysis.core_pan.compare_to_core_pan(model,pcmodel,based_on)
   report.save(dir)


@analyse.command()
@click.argument('modelpaths', nargs=-1, type=click.Path(exists=True))
# @TODO 
#    extend / add different comparison options
@click.option('--type','-t',required=False,type=click.Choice(['sboterm']), multiple=True,
              show_default=True,  default=[], help='Type of comparison to be performed.')
@click.option('--all',required=False, is_flag=True, show_default=True, default=False, 
              help='Shortcut to run all comparisons. Overwrites input of type.')
@click.option('-d', '--dir', required=False, type=click.Path(file_okay=False), show_default=True, default='comparison', help='Path to the output dir.')
@click.option('-c','--colour', required=False, type=str, default='Paired', show_default=True, help='Name of Matplotlib colour palette for the plot.')
def compare(modelpaths,type,all,dir,colour):
   """Compare models.
   """
   Path(dir).mkdir(parents=True, exist_ok=True)
   # compare SBOterms
   if all or 'sboterm' in type:
      models = rg.utility.io.load_model(list(modelpaths),'libsbml')
      fig = rg.analysis.comparison.plot_rea_sbo_multiple(models,color_palette=colour)
      fig.savefig(Path(dir, 'sboterm'), dpi=400, bbox_inches='tight')


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
@click.option('-n', '--namespace', required=False, type=click.Choice(['BiGG']), show_default=True, default='BiGG', help='Namespace to use for the model.')
@click.option('-d', '--dir', required=False, type=click.Path(), show_default=True, default='', help='Path to the output dir.')
@click.option('-c', '--colors', required=False, type=str, show_default=True, default='YlGn', help='Abbreviation of a matplotlib colour palette.')
def simulate(modelpaths,media,namespace,dir,colors):
   """Simulate the growth of the given model vs. media.
   """
   report = rg.analysis.growth.growth_analysis(modelpaths, media, namespace,retrieve='report')
   report.save(to=dir, color_palette=colors)


# @TEST
@growth.command()
@click.argument('modelpath', type=click.Path(exists=True))
@click.option('-m','--media',required=True,type=click.Path(exists=True), help='Path to a media config file.')
@click.option('-n', '--namespace', required=False, type=click.Choice(['BiGG']), show_default=True, default='BiGG', help='Namespace to use for the model.')
@click.option('-d', '--dir', required=False, type=click.Path(), show_default=True, default='', help='Path to the output dir.')
@click.option('-c', '--colors', required=False, type=str, show_default=True, default='YlGn', help='Abbreviation of a matplotlib colour palette.')
def auxotrophies(modelpath,media,namespace,dir,colors):
   """Test for auxotrophies for the 20 proteinogenic amino acids.
   """
   medialist, supps = rg.analysis.growth.read_media_config(media)
   model = rg.utility.io.load_model(modelpath, 'cobra')
   report = rg.analysis.growth.test_auxotrophies(model,medialist,namespace)
   report.save(dir,colors)

#@TODO --substances
@growth.command()
@click.argument('modelpath', type=click.Path(exists=True))
@click.option('-e','--element', type=str, required=True, help='Element to perform the source test with. Needs to be a valid chemical, elemental Symbol.')
@click.option('-s','--substances', type=str, required=False, show_default=True, default=None, multiple=True, help='Add substances from the database to test against. If none are given, resrs against the while database.')
@click.option('-m', '--medium', required=False, type=str, show_default=True, default=None, help='Name of a medium from the database to use for the testing. If not given, uses the one from the model.')
@click.option('-n', '--namespace', required=False, type=click.Choice(['BiGG']), show_default=True, default='BiGG', help='Namespace to use for the model.')
@click.option('-d', '--dir', required=False, type=click.Path(), default='', help='Path to the output dir.')
@click.option('-c', '--colors', required=False, type=str, show_default=True, default='YlGn', help='Abbreviation of a matplotlib colour palette.')
def sources(modelpath, element, substances, medium, namespace, dir, colors):
   """Test growth on different sources for an element.
   """ 
   model = rg.utility.io.load_model(modelpath,'cobra')
   report = rg.analysis.growth.test_growth_with_source(model,element,substances,medium,namespace)
   report.save(dir,color_palette=colors)


@growth.command()
@click.argument('modelpath', type=click.Path(exists=True))
@click.option('-o','--objective', required=False, type=click.Choice(['flux','medium','exchanges']), show_default=True, default='flux',help='Set the type of minimal medium to be calculated: minimal fluxes, minimal compounds or based on the model exchange reactions.')
@click.option('-r', '--growth-rate', required=False, type=float, show_default=True, default=0.5, help='Minimal rate to be reached by the minimal medium. The smaller the rate the more expensive the calculation.')
@click.option('-d', '--dir', required=False, type=click.Path(), show_default=True, default='', help='Path to the output dir.')
def minimal_medium(modelpath,objective, growth_rate, dir):
   """Calculate the minimal medium of a model.

   Either set the fluxes minimal, the number of compounds based on current medium or the number of exchange reactions.
   """
    
   model = rg.utility.io.load_model(modelpath,'cobra')
   medium = rg.analysis.growth.model_minimal_medium(model, objective, growth_rate)
   medium.export_to_file(dir=dir)
