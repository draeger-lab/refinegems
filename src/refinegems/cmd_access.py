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
# @TODO gaps group still for the old gapfill - rewrite or delete
# @TEST help is displayed alright but untested
@cli.group()
def gaps():
   """Find and fill gaps in a model."""

@gaps.command(show_constraints=True)
@cloup.argument('alg', type=click.Choice(['KEGG','BioCyc','Gene']),
                help='Type of automated gap-gapfilling algorithm, that shall be used.')
@cloup.argument('modelpath', type=click.Path(exists=True, dir_okay=False), 
                help='Path to the model.')
@cloup.option_group(
   "General options",
   cloup.option('-o','--outdir', type=click.Path(exists=True, file_okay=False), 
                default='./', show_default=True, 
                help='Path to a directory to write the output to.'),
   cloup.option('-f','--fill', is_flag=True, type=bool, default=False,
                help='If True, tries to fill the gaps in the model.'),
   cloup.option('--fc','--formula-check', type=click.Choice(['none','existence','wildcard','strict']),
                default='existence', show_choices=True, show_default=True,
                help='Set the filter for which metabolite formulas are valid to be added to the model.'),
   cloup.option('--no-dna', is_flag=True, type=bool, default=False, 
                help='Exclude DNA reactions (name-based) from being added to the model.'),
   cloup.option('--no-rna', is_flag=True, type=bool, default=False, 
                help='Exclude RNA reactions (name-based) from being added to the model.'),
   cloup.option('-p','--idprefix', type=str, default='refineGEMs', show_default=True,
                help='Prefix for the random IDs, if an ID does not exists for the given namespace.'),
   cloup.option('-n','--namespace', type=click.Choice(['BiGG']), default='BiGG',
                show_default=True, help='Namespace used in the model.')
)
@cloup.option_group(
   "KEGG required parameters",
   "Parameters required when running the KEGG gap-gapfilling algorithmn",
   cloup.option('--orgid', type=str, help='KEGG organism ID'),
   constraint = cloup.constraints.If(cloup.constraints.Equal('alg','KEGG'), 
                                     then=cloup.constraints.require_all)
)
@cloup.option_group(
   "BioCyc required parameters",
   "Parameters required when running the KEGG gap-gapfilling algorithmn",
   cloup.option('--gt','--genetable', type=click.Path(exists=True, dir_okay=False),
                help='Path to the BioCyc gene smart table.'),
   cloup.option('--rt','--reactable', type=click.Path(exists=True, dir_okay=False),
                help='Path to the BioCyc gene smart table.'),
   cloup.option('--gff-bc', type=click.Path(exists=True, dir_okay=False),
                help='Path to the GFF.'),
   constraint = cloup.constraints.If(cloup.constraints.Equal('alg','BioCyc'), 
                                     then=cloup.constraints.require_all)
)
@cloup.option_group(
   "Gene required parameters",
   "Parameters required when running the GeneGapFiller algorithm",
   cloup.option('--gff-g', type=click.Path(exists=True, dir_okay=False),
                help='Path to the GFF.'),
   constraint = cloup.constraints.If(cloup.constraints.Equal('alg','Gene'), 
                                     then=cloup.constraints.require_all)
)
@cloup.option_group(
   "Gene optional parameters",
   "Optional / conditionally interdependant parameters for the gene gap-gapfilling algorithm",
   cloup.option('--prot-prefix', type=str, default='refineGEMs',
                show_default = True,
                help='Prefix for pseudo-protein IDs.'),
   cloup.option('--mail', type=str, default=None, help='Mail address for NCBI requests.'),
   cloup.option('--ncbi','--check-ncbi', is_flag=True, default=False,
                help='Enable searching protein IDs in NCBI. This increases the runtime significantly.'),
   cloup.option('--fasta', type=click.Path(exists=True, dir_okay=False), default=None,
                help='Path to the protein FASTA of the model.'),
   cloup.option('--dmnd-db', type=click.Path(exists=True, dir_okay=False), default=None,
                 help = 'Path to the SwissProt DIAMOND database.'),
   cloup.option('--sp-map','--swissprot-mapping', type=click.Path(exists=True, dir_okay=False),
                default=None, help='Path to the SwissProt mapping file (ID against EC and BRENDA)'),
   cloup.option('-s','--sensitivity', type=click.Choice(['sensitive', 'more-sensitive', 'very-sensitive','ultra-sensitive']),
                default='more-sensitive', show_default=True, show_choices=True,
                help='Sensitivity mode for running DIAMOND.'),
   cloup.option('--cov', type=float, default=90.0, show_default=True,
                help='Coverage value (passed to DIAMOND)'),
   cloup.option('--pid', type=float, default=95.0, show_default=True,
                help='Percentage identity threshold value for filtering DIAMOND results.'),
   cloup.option('-t','--threads',type=int, default=2, show_default=True, 
                help='Number of threads to be used by DIAMOND.'), 
)
@cloup.constraints.constraint(cloup.constraints.If(cloup.constraints.IsSet('ncbi'), then=cloup.constraints.require_all), ['mail'])
@cloup.constraints.constraint(cloup.constraints.If(cloup.constraints.AnySet('fasta','dmnd_db','sp_map'),then=cloup.constraints.require_all), ['fasta','dmnd_db','sp_map'])
@cloup.constraints.constraint(cloup.constraints.If(cloup.constraints.AnySet('sensitivity','cov','pid','threads'),then=cloup.constraints.require_all), ['fasta','dmnd_db','sp_map'])
def automated_gapfill(alg,modelpath,outdir,fill,
                      formula_check, no_dna, no_rna,
                      idprefix, namespace,
                      orgid,
                      genetable,reactable,gff_bc,
                      gff_g,
                      prot_prefix, mail, ncbi, 
                      fasta, dmnd_db, sp_map, sensitivity, cov, pid, threads):
   
   # @TODO incomplete 
   pass 
   
   cmodel = rg.utility.io.load_model(modelpath, 'cobra')
   model = rg.utility.io.load_model(modelpath, 'libsbml')
   
   # set class instance
   match alg:
      case 'KEGG':
         gapfiller = rg.classes.gapfill.KEGGapFiller(orgid)
         # find gaps
         gapfiller.find_missing_genes(model)
         gapfiller.find_missing_reactions(cmodel)
      case 'BioCyc':
         gapfiller = rg.classes.gapfill.BioCycGapFiller(genetable,
                                                        reactable,
                                                        gff_bc)
         # find gaps
         gapfiller.find_missing_genes(model)
         gapfiller.find_missing_reactions(cmodel)
      case 'Gene':
         gapfiller = rg.classes.gapfill.GeneGapFiller()
         # find gaps
         gapfiller.find_missing_genes(gff_g,model)
         gapfiller.find_missing_reactions(cmodel, prot_prefix, mail, ncbi, fasta, dmnd_db, sp_map,
                                 sensitivity, cov, pid, threads)
      case _:
         mes = f'Unknown option for algorthmn type: {alg}'
         raise ValueError(mes)
   # fill gaps
   if fill:
      model = gapfiller.fill_model(model, 
                                   formula_check=formula_check,
                                   exclude_dnae=no_dna, exclude_rna=no_rna,
                                   idprefix=idprefix, namespace=namespace)
      # save model
      write_model_to_file(model, Path(outdir, 'gapfilled_model.xml'))
      
   # @TODO report stats
   # @TODO report manual curation
   
   


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
