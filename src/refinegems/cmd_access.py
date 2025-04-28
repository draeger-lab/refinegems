"""Entry points to the code from the command line."""

__author__ = "Carolin Brune, Gwendolyn O. Döbel"

################################################################################
# Requirements
################################################################################

from pathlib import Path

import click
import cloup
import pandas as pd

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
    """Set-up tools, folder structure and more for running the program."""


# download important stuff
# ------------------------
@setup.command()
@click.option(
    "--filename",
    "-f",
    default="config.yaml",
    type=str,
    show_default=True,
    help="Name (path) to save the config file under.",
)
@click.option(
    "--type",
    "-t",
    show_default=True,
    default="media",
    type=click.Choice(["media", "refinegems"]),
    help="Type of config file to download. Either a file for media configuration or a file to run a refinement pipeline.",
)
def config(filename, type):
    """Download a configuration file (.yaml).

    Download a configuration file to edit for running a specific function
    or provide the media configuration file.
    """
    rg.utility.set_up.download_config(filename, type)


@setup.command()
@click.argument("downloadtype", type=click.Choice(["SwissProt_gapfill"]))
@click.option(
    "-d",
    "--dir",
    required=False,
    type=click.Path(),
    show_default=True,
    default="",
    help="Path to the output dir.",
)
@click.option(
    "-c",
    "--chunksize",
    required=False,
    type=int,
    show_default=True,
    default=10,
    help="Chunksizre for the download in kB. Defaults to 10.",
)
def data(downloadtype, dir, chunksize):
    """Download file(s) needed for a given functionality of the toolbox."""
    dltype = downloadtype.replace("_", " ")
    rg.utility.set_up.download_url(dltype, dir, chunksize)


# Get ID mapping for GeneProducts
# -------------------------------
@setup.command()
@click.argument("modelpath", type=str)
@click.option(
    "-g",
    "--gff-paths",
    required=False,
    show_default=True,
    default=None,
    type=click.Path(exists=True, dir_okay=False),
    multiple=True,
    help="Path(s) to GFF file(s). Allowed GFF formats are: RefSeq, NCBI and Prokka. ",
)
@click.option(
    "-e",
    "--email",
    required=False,
    show_default=True,
    default=None,
    type=str,
    help="E-mail for NCBI queries. This is only used when --mapping-tbl-file is not provided.",
)
@click.option(
    "-t",
    "--lt",
    "--contains-locus-tags",
    show_default=True,
    default=False,
    type=bool,
    help="Specifies if provided model has locus tags within the label tag if set to True.",
)
@click.option(
    "-o",
    "--outdir",
    required=False,
    type=click.Path(exists=False, dir_okay=True, file_okay=False, writable=True, path_type=Path),
    default=Path(""),
    show_default=True,
    help="Path to a directory to write the output files to.",
)
def geneproduct_mapping_table(modelpath, gff_paths, email, lt, outdir):
    """Generates ID mapping file for GeneProducts"""
    model = rg.utility.io.load_model(modelpath, "libsbml")
    rg.utility.entities.get_gpid_mapping(
        model, gff_paths, email, lt, outdir
    )


# build a pan-core model
# ----------------------
@setup.command()
@click.argument("models", nargs=-1, type=click.Path())
@click.option(
    "-o",
    "--based-on",
    required=False,
    type=click.Choice(["id"]),
    show_default=True,
    default="id",
    help="Option on how to combine the models.",
)
@click.option(
    "-n",
    "--name",
    required=False,
    type=str,
    show_default=True,
    default="pan-core-model",
    help="Name of the new pan-core model.",
)
@click.option(
    "-g",
    "--keep-genes",
    is_flag=True,
    show_default=True,
    default=False,
    help="If set, the genes are kept in the pan-core model, otherwise they are deleted.",
)
@click.option(
    "--rcomp",
    "--resolve-compartments",
    is_flag=True,
    help="If set, tries to standardise the compartment names to the c,p,e,... namespace.",
)
@click.option(
    "-d",
    "--dir",
    required=False,
    type=click.Path(),
    show_default=True,
    default="",
    help="Path to the output dir.",
)
def build_pancore(models, based_on, name, keep_genes, rcomp, dir):
    """Build a pan-core model."""
    models = list(models)
    # construc the core-pan model
    pancore_mod = rg.analysis.core_pan.generate_core_pan_model(
        models, based_on, name, not keep_genes
    )
    # resolve compartment names
    if rcomp:
        rg.analysis.core_pan.resolve_compartment_names(pancore_mod)
    # save the model
    rg.utility.io.write_model_to_file(pancore_mod, str(Path(dir, name + ".xml")))


# --------------------
# work on the database
# --------------------


@cli.group()
def database():
    """Access, curate, etc. the in-build database."""


@database.command()
def initialise():
    """Initialise or update the database."""
    rg.utility.databases.initialise_database()


@database.command()
@click.argument("databasename", type=click.Choice(["MetaNetX"]))
@click.option(
    "-c",
    "--chunksize",
    show_default=True,
    default=1,
    help="Chunksize for the download in kB. Defaults to 1.",
)
def add_namespace(databasename, chunksize):
    """Add or update table for a given namespace/database to the database."""
    match databasename:
        case "MetaNetX":
            rg.utility.databases.update_mnx_namespaces(chunksize=chunksize)
        case _:
            mes = f"Unknown database name: {databasename}"
            raise ValueError(mes)


@database.command()
def reset():
    """Reset the database by removing additionally downloaded tables."""
    rg.utility.databases.reset_database()


# -------------------
# all about the media
# -------------------


@cli.group()
def media():
    """Access the media part of the database."""


# Handle medium / media database
# ------------------------------
@media.command()
@click.option(
    "--list",
    is_flag=True,
    show_default=True,
    default=False,
    help="List the names of the media in the database.",
)
# @click.option('--copy', is_flag=False, flag_value='./media.csv', default=None, type=click.Path(exists=False),
# help='Produce and save a copy of the media database.')
def info(list):  # ,copy
    """Access information about the in-build media database.

    Can be used to check the available media
    """
    if list:
        possible_media = rg.utility.io.load_a_table_from_database("medium", False)[
            "name"
        ].to_list()
        media = "|".join(possible_media)
        print(media)


"""
    if copy:
        db = specimen.classes.medium.load_media_db()
        specimen.classes.medium.save_db(db,click.format_filename(copy))
"""


# --------------
# Curate a model
# --------------
@cli.group()
def curate():
    """Curate a model"""


# Polish annotations
# ------------------
@curate.command()
@click.argument("modelpath", type=str)
@click.option(
    "-n",
    "--new-pattern",
    required=False,
    show_default=True,
    default=True,
    type=bool,
    help="Specifies if new pattern `database_prefix:local_identifier` for CURIEs should be used if set to True. ",
)
@click.option(
    "-o",
    "--outdir",
    required=False,
    type=click.Path(exists=False, dir_okay=True, file_okay=False, writable=True, path_type=Path),
    default=Path(""),
    show_default=True,
    help="Path to a directory to write the output files to.",
)
def annotations(modelpath, new_pattern, outdir):
    """Polishes annotations

    Changes qualifiers and annotations to be MIRIAM-compliant
    """
    model = rg.utility.io.load_model(modelpath, "libsbml")
    model = rg.curation.curate.polish_annotations(model, new_pattern, outdir)
    rg.utility.io.write_model_to_file(
        model, str(Path(outdir, f"{model.getId()}_after_annots_polishing.xml"))
    )


# Perform all steps to polish a model
# -----------------------------------
@curate.command()
@cloup.argument(
    "modelpath", type=click.Path(exists=True, dir_okay=False), help="Path to the model."
)
@cloup.option_group(
    "General options",
    cloup.option(
        "-i",
        "--id-db",
        required=False,
        show_default=True,
        default="BIGG",
        type=str,
        help="Main database where identifiers in model come from",
    ),
    cloup.option(
        "-m",
        "--mtf",
        "--mapping-tbl-file",
        required=False,
        show_default=True,
        default=None,
        type=click.Path(exists=True, dir_okay=False),
        help="Path to a file containing a mapping table with columns `model_id | X...` where `X` can be `REFSEQ`, `NCBI`, `locus_tag` or `UNCLASSIFIED`. The table can contain all of the ``X`` columns or at least one of them. ",
    ),
    cloup.option(
        "-l",
        "--lab-strain",
        required=False,
        show_default=True,
        default=False,
        type=bool,
        help="Specifies if a strain from no database was provided and thus has only homolog mappings, if set to True. ",
    ),
    cloup.option(
        "-k",
        "--kid",
        "--kegg-organism-id",
        required=False,
        show_default=True,
        default=None,
        type=str,
        help="KEGG organism identifier if available.",
    ),
    cloup.option(
        "-r",
        "--rd",
        "--reaction-direction",
        required=False,
        type=click.Path(exists=True, dir_okay=False),
        default=None,
        show_default=True,
        help="Path to a CSV file containing the BioCyc smart table with the columns `Reactions (MetaCyc ID) | EC-Number | KEGG reaction | METANETX | Reaction-Direction`.",
    ),
    cloup.option(
        "-o",
        "--outdir",
        required=False,
        type=click.Path(exists=False, dir_okay=True, file_okay=False, writable=True, path_type=Path),
        default=Path(""),
        show_default=True,
        help="Path to a directory to write the output files to.",
    ),
)
@cloup.option_group(
    "Parameters required when --mapping-tbl-file is not provided",
    cloup.option(
        "-g",
        "--gff-paths",
        required=False,
        show_default=True,
        default=None,
        type=click.Path(exists=True, dir_okay=False),
        multiple=True,
        help="Path(s) to GFF file(s). Allowed GFF formats are: RefSeq, NCBI and Prokka.",
    ),
    cloup.option(
        "-e",
        "--email",
        required=False,
        show_default=True,
        default=None,
        type=str,
        help="E-mail for NCBI queries.",
    ),
    cloup.option(
        "-t",
        "--lt",
        "--contains-locus-tags",
        show_default=True,
        default=False,
        type=bool,
        help="Specifies if provided model has locus tags within the label tag if set to True.",
    ),
)
def model(
    modelpath,
    id_db,
    mtf,
    gff_paths,
    email,
    lt,
    lab_strain,
    kid,
    rd,
    outdir,
):
    """Completes all steps to polish a model created by an automatic reconstruction pipeline.

    Extends annotatations, cleans up the notes fields and changes qualifiers and annotations to be MIRIAM-compliant.

    (Tested for models having either BiGG or VMH identifiers.)
    """
    model = rg.utility.io.load_model(modelpath, "libsbml")
    model = rg.curation.curate.polish_model(
        model,
        id_db,
        mtf,
        gff_paths,
        email,
        lt,
        lab_strain,
        kid,
        rd,
        outdir,
    )
    rg.utility.io.write_model_to_file(
        model, str(Path(outdir, f"{model.getId()}_after_polishing.xml"))
    )

# --------------
# Refine a model
# --------------
@cli.group()
def refine():
    """Refine a model. Includes steps like biomass, charges, SBO annotation, reaction direction correction and addition
    of Pathways and further gene product annotations.
    """


# Find and fill gaps in a model automatically/Fill gaps with manually created tables
# ----------------------------------------------------------------------------------

@refine.command(show_constraints=True)
@cloup.argument(
    "alg",
    type=click.Choice(["KEGG", "BioCyc", "Gene"]),
    help="Type of automated gapfilling algorithm, that shall be used.",
)
@cloup.argument(
    "modelpath", type=click.Path(exists=True, dir_okay=False), help="Path to the model."
)
@cloup.option_group(
    "General options",
    cloup.option(
        "-o",
        "--outdir",
        type=click.Path(exists=False, dir_okay=True, file_okay=False, writable=True, path_type=Path),
        default="./",
        show_default=True,
        help="Path to a directory to write the output to.",
    ),
    cloup.option(
        "-f",
        "--fill",
        is_flag=True,
        type=bool,
        default=False,
        help="If True, tries to fill the gaps in the model.",
    ),
    cloup.option(
        "-fc",
        "--formula-check",
        type=click.Choice(["none", "existence", "wildcard", "strict"]),
        default="existence",
        show_choices=True,
        show_default=True,
        help="Set the filter for which metabolite formulas are valid to be added to the model.",
    ),
    cloup.option(
        "--no-dna",
        is_flag=True,
        type=bool,
        default=False,
        help="Exclude DNA reactions (name-based) from being added to the model.",
    ),
    cloup.option(
        "--no-rna",
        is_flag=True,
        type=bool,
        default=False,
        help="Exclude RNA reactions (name-based) from being added to the model.",
    ),
    cloup.option(
        "-p",
        "--idprefix",
        type=str,
        default="refineGEMs",
        show_default=True,
        help="Prefix for the random IDs, if an ID does not exists for the given namespace.",
    ),
    cloup.option(
        "-n",
        "--namespace",
        type=click.Choice(["BiGG"]), 
        default="BiGG",
        show_default=True,
        help="Namespace used in the model.",
    ),
    cloup.option(
        "-r",
        "--report",
        is_flag=True,
        type=bool,
        default=False,
        help="If True, saves statistics and genes/reactions for manual curation.",
    ),
)
@cloup.option_group(
    "KEGG required parameters",
    "Parameters required when running the KEGG gapfilling algorithmn",
    cloup.option("--orgid", type=str, help="KEGG organism ID"),
    constraint=cloup.constraints.If(
        cloup.constraints.Equal("alg", "KEGG"), then=cloup.constraints.require_all
    ),
)
@cloup.option_group(
    "BioCyc required parameters",
    "Parameters required when running the BioCyc gapfilling algorithmn",
    cloup.option(
        "-gt",
        "--genetable",
        type=click.Path(exists=True, dir_okay=False),
        help="Path to the BioCyc gene smart table.",
    ),
    cloup.option(
        "-rt",
        "--reactable",
        type=click.Path(exists=True, dir_okay=False),
        help="Path to the BioCyc reaction smart table.",
    ),
    cloup.option(
        "--gff-bc",
        type=click.Path(exists=True, dir_okay=False),
        help="Path to the GFF.",
    ),
    constraint=cloup.constraints.If(
        cloup.constraints.Equal("alg", "BioCyc"), then=cloup.constraints.require_all
    ),
)
@cloup.option_group(
    "Gene required parameters",
    "Parameters required when running the GeneGapFiller algorithm",
    cloup.option(
        "--gff-g", type=click.Path(exists=True, dir_okay=False), help="Path to the GFF."
    ),
    constraint=cloup.constraints.If(
        cloup.constraints.Equal("alg", "Gene"), then=cloup.constraints.require_all
    ),
)
@cloup.option_group(
    "Gene optional parameters",
    "Optional / conditionally interdependant parameters for the gene gap-gapfilling algorithm",
    cloup.option(
        "--prot-prefix",
        type=str,
        default="refineGEMs",
        show_default=True,
        help="Prefix for pseudo-protein IDs.",
    ),
    cloup.option(
        "--mail", type=str, default=None, help="Mail address for NCBI requests.",
    ),
    cloup.option(
        "-ncbi",
        "--check-ncbi",
        is_flag=True,
        default=False,
        help="Enable searching protein IDs in NCBI. This increases the runtime significantly.",
    ),
    cloup.option(
        "--db-type",
        type = click.Choice(['swissprot','user']),
        default='swissprot',
        show_default=True,
        show_choices=True,
        help="Database to search against. Choose 'swissprot' for SwissProt or 'user' for a user defined database.",
    ),  
    cloup.option(
        "--fasta",
        type=click.Path(exists=True, dir_okay=False),
        default=None,
        help="Path to the protein FASTA of the model.",
    ),
    cloup.option(
        "--dmnd-db",
        type=click.Path(exists=True, dir_okay=False),
        default=None,
        help="Path to the SwissProt DIAMOND database.",
    ),
    cloup.option(
        "-db-map",
        "--db-mapping",
        type=click.Path(exists=True, dir_okay=False),
        default=None,
        help="Path to the SwissProt or User defined mapping file (SwissProt: ID against EC and BRENDA), User: ID against EC",
    ),
    cloup.option(
        "-s",
        "--sensitivity",
        type=click.Choice(
            ["sensitive", "more-sensitive", "very-sensitive", "ultra-sensitive"]
        ),
        default="more-sensitive",
        show_default=True,
        show_choices=True,
        help="Sensitivity mode for running DIAMOND.",
    ),
    cloup.option(
        "--cov",
        type=float,
        default=90.0,
        show_default=True,
        help="Coverage value (passed to DIAMOND)",
    ),
    cloup.option(
        "--pid",
        type=float,
        default=95.0,
        show_default=True,
        help="Percentage identity threshold value for filtering DIAMOND results.",
    ),
    cloup.option(
        "--threshold-add-reacs",
        type=int,
        default=5,
        show_default=True,
        help="Maximal number of reactions for one EC number mapping for it to be considered successful and to be added to the model.",
    ),
    cloup.option(
        "-t",
        "--threads",
        type=int,
        default=2,
        show_default=True,
        help="Number of threads to be used by DIAMOND.",
    ),
)
@cloup.constraints.constraint(
    cloup.constraints.If(
        cloup.constraints.IsSet("check_ncbi"), then=cloup.constraints.require_all
    ),
    ["mail"],
)
@cloup.constraints.constraint(
    cloup.constraints.If(
        cloup.constraints.AnySet("db_type", "fasta", "dmnd_db", "db_mapping"),
        then=cloup.constraints.require_all,
    ),
    ["db_type","fasta", "dmnd_db", "db_mapping"],
)
# @cloup.constraints.constraint(cloup.constraints.If(cloup.constraints.AnySet('sensitivity','cov','pid','threads'),then=cloup.constraints.require_all), ['fasta','dmnd_db','sp_map'])
def automated_gapfill(
    alg,
    modelpath,
    outdir,
    fill,
    formula_check,
    no_dna,
    no_rna,
    idprefix,
    namespace,
    report,
    orgid,
    genetable,
    reactable,
    gff_bc,
    gff_g,
    prot_prefix,
    mail,
    check_ncbi,
    db_type,
    fasta,
    dmnd_db,
    db_mapping,
    sensitivity,
    cov,
    pid,
    threshold_add_reacs,
    threads,
):

    cmodel = rg.utility.io.load_model(modelpath, "cobra")
    model = rg.utility.io.load_model(modelpath, "libsbml")

    # set class instance
    match alg:
        case "KEGG":
            gapfiller = rg.classes.gapfill.KEGGapFiller(orgid)
            # find gaps
            gapfiller.find_missing_genes(model)
            gapfiller.find_missing_reactions(cmodel, threshold_add_reacs)
        case "BioCyc":
            gapfiller = rg.classes.gapfill.BioCycGapFiller(genetable, reactable, gff_bc)
            # find gaps
            gapfiller.find_missing_genes(model)
            gapfiller.find_missing_reactions(cmodel)
        case "Gene":
            gapfiller = rg.classes.gapfill.GeneGapFiller()
            # find gaps
            gapfiller.find_missing_genes(gff_g, model)
            kwargs = {
                "outdir": outdir,
                "sens": sensitivity,
                "cov": cov,
                "t": threads,
                "pid": pid,
            }
            gapfiller.find_missing_reactions(
                model = cmodel,
                prefix = prot_prefix,
                type_db = db_type,
                fasta = fasta,
                dmnd_db = dmnd_db,
                map_db = db_mapping,
                mail = mail,
                check_NCBI = check_ncbi,
                threshold_add_reacs = threshold_add_reacs,
                **kwargs,
            )
        case _:
            mes = f"Unknown option for algorithm type: {alg}"
            raise ValueError(mes)
    # fill gaps
    if fill:
        model = gapfiller.fill_model(
            model,
            formula_check=formula_check,
            exclude_dna=no_dna,
            exclude_rna=no_rna,
            idprefix=idprefix,
            namespace=namespace,
        )
        # save model
        write_model_to_file(model, str(Path(outdir, "gapfilled_model.xml")))
    # report statistics
    if report:
        gapfiller.report(outdir)



@refine.command()
@click.argument("modelpath", type=click.Path(exists=True))
@click.option(
    "-c",
    "--cycles",
    show_default=True,
    default=10,
    type=int,
    help="Maximal number of optimisation cycles that will be run.",
)
@click.option(
    "--outfile",
    "-o",
    required=False,
    type=click.Path(),
    show_default=True,
    default="biomass_corrected.xml",
    help="Path to save the corrected model under.",
)
def biomass(modelpath, cycles, outfile):
    """Changes the biomass reaction to be consistent"""
    model = rg.utility.io.load_model(modelpath, "cobra")
    corrected = rg.curation.biomass.check_normalise_biomass(model, cycles)
    if corrected:
        rg.utility.io.write_model_to_file(model, outfile)


@refine.command()
@click.argument("modelpath", type=click.Path(exists=True))
@click.option(
    "--dir",
    "-d",
    required=False,
    type=click.Path(),
    show_default=True,
    default="",
    help="Directory to save the output to.",
)
def charges(modelpath, dir):
    """Changes the charges in a model by comparison to the ModelSEED database"""
    model = rg.utility.io.load_model(modelpath, "libsbml")
    corr, charg = rg.curation.charges.correct_charges_modelseed(model)
    rg.utility.io.write_model_to_file(corr, str(Path(dir, "corrected_model.xml")))
    charg_tab = pd.DataFrame.from_dict({k: pd.Series(v) for k, v in charg.items()})
    charg_tab.to_csv(Path(dir, "multiple_charges.csv"), sep=";")


@refine.command()
@click.argument("modelpath", type=click.Path(exists=True))
@click.argument("data", type=click.Path(exists=True))
@click.option(
    "--dir",
    "-d",
    required=False,
    type=click.Path(),
    show_default=True,
    default="",
    help="Directory to save the output to.",
)
def direction(modelpath, data, dir):
    """Checks & if necessary, corrects the direction for each reaction in a model"""
    model = rg.utility.io.load_model(modelpath, "cobra")
    corr = rg.curation.polish.check_direction(model, data)
    rg.utility.io.write_model_to_file(
        corr, str(Path(dir, f"{model.id}_corrected_model.xml"))
    )


# egcs
@refine.command()
@click.argument("modelpath", type=click.Path(exists=True))
@click.option(
    "--solver",
    "-s",
    required=False,
    type=click.Choice(["greedy"]),
    show_default=True,
    default=None,
    multiple=False,
    help="Type of solver for the EGCs.",
)
@click.option(
    "--namespace",
    "-n",
    required=False,
    type=click.Choice(["BiGG"]),
    show_default=True,
    default="BiGG",
    multiple=False,
    help="Namespace of the model.",
)
@click.option(
    "--compartment",
    "-c",
    required=False,
    type=str,
    show_default=True,
    default="c,e",
    help="Compartments to check, separated by comma only.",
)
@click.option(
    "--outfile",
    "-o",
    required=False,
    type=click.Path(),
    show_default=True,
    default="fixed_egcs.xml",
    help="Path to save edited model to, if solver has been set.",
)
def egcs(modelpath, solver, namespace, compartment, outfile):
    """Identify and optionally solve EGCs."""
    model = rg.utility.io.load_model(modelpath, "cobra")
    compartment = compartment.split(",")
    match solver:
        case "greedy":
            solver = rg.classes.egcs.GreedyEGCSolver()
            solution = solver.solve_egcs(model, namespace, compartment)
            if solution:
                for k, v in solution.items():
                    print(k + ": " + v)
                rg.utility.io.write_model_to_file(model, outfile)
        case _:
            solver = rg.classes.egcs.EGCSolver()
            egcs, vals = solver.find_egcs(model, True, namespace, compartment)
            print(f"EGCs: \n{str(egcs)}")
            print(f"Objective values: \n{str(vals)}")


# Annotation-related clean-up
# & Additional annotations
# ---------------------------
@refine.group()
def annot():
    """Add annotations to your model."""


@annot.command()
@click.argument("modelpath", type=click.Path(exists=True))
@click.option(
    "--dir",
    "-d",
    required=False,
    type=click.Path(),
    show_default=True,
    default="",
    help="Directory to save the output to.",
)
def sboterms(modelpath, dir):
    """Calls SBOannotator to enhance the SBO terms in a model"""
    model = rg.utility.io.load_model(modelpath, "libsbml")
    SBOannotated = rg.utility.connections.run_SBOannotator(model)
    rg.utility.io.write_model_to_file(
        SBOannotated, str(Path(dir, "SBOannotated_model.xml"))
    )


@annot.command()
@click.argument("modelpath", type=click.Path(exists=True))
@click.option(
    "-d", "--dir", required=False, default="", help="Path to the output directory"
)
def pathways(modelpath, dir):
    """Add KEGG pathways as groups to a model"""
    model, missing = rg.curation.pathways.kegg_pathways(modelpath)
    with open(Path(dir, "reac_wo_kegg_pathway_groups.txt"), "w") as outfile:
        # save reactions with missing groups
        for line in missing:
            outfile.write(f"{line}\n")
    # save model
    write_model_to_file(model, str(Path(dir, "model_with_added_KeggPathwayGroups.xml")))


# -----------------------------------------------
# Analyse a model/multiple models including a comparison
# -----------------------------------------------
@cli.group()
def analyse():
    """Analyse a model by testing for auxotrophies and growth on different media along with finding EGCs and looking at
    overall model statistics.
    """


@analyse.command()
@click.argument("modelpath", type=click.Path(exists=True))
@click.option(
    "-s",
    "--score-only",
    is_flag=True,
    show_default=True,
    default=False,
    help="Specifies if memote is only run to return the score",
)
@click.option(
    "-f",
    "--file",
    required=False,
    type=str,
    show_default=True,
    default="memote.html",
    help="Name or path to save the output to. Only relevant if -s is not set.",
)
def memote(modelpath, score_only, file):
    """Perform a memote analysis."""
    model = rg.utility.io.load_model(modelpath, "cobra")
    if score_only:
        score = rg.utility.connections.get_memote_score(
            rg.utility.connections.run_memote(model, type="json", return_res=True)
        )
        print(f"The Memote score is {score}.")
    else:
        rg.utility.connections.run_memote(model, save_res=file)


@analyse.command()
@click.argument("modelpath", type=click.Path(exists=True))
@click.option(
    "-d",
    "--dir",
    required=False,
    type=click.Path(),
    show_default=True,
    default="",
    help="Path to the output dir.",
)
@click.option(
    "-c",
    "--colors",
    required=False,
    type=str,
    show_default=True,
    default="YlGn",
    help="Abbreviation of a matplotlib colour palette.",
)
def pathways(modelpath, dir, colors):
    """Analysis of pathways contained in a model"""
    model = rg.utility.io.load_model(modelpath, "cobra")
    report = rg.curation.pathways.kegg_pathway_analysis(model)
    report.save(dir, colors)


@analyse.command()
@click.argument("modelpath", type=click.Path(exists=True))
@click.option(
    "-d",
    "--dir",
    required=False,
    type=click.Path(),
    show_default=True,
    default="",
    help="Path to the output dir.",
)
@click.option(
    "-c",
    "--colors",
    required=False,
    type=str,
    show_default=True,
    default="YlGn",
    help="Abbreviation of a matplotlib colour palette.",
)
def stats(modelpath, dir, colors):
    """Generate a report on the statistics of a model."""
    model = rg.utility.io.load_model(modelpath, "cobra")
    report = rg.classes.reports.ModelInfoReport(model)
    report.save(dir, colors)


@analyse.command()
@click.argument("modelpath", type=click.Path(exists=True))
@click.argument("pcpath", type=click.Path(exists=True))
@click.option(
    "-b",
    "--based-on",
    type=click.Choice(["id"]),
    required=False,
    show_default=True,
    default="id",
    help="Option on how to compare the models.",
)
@click.option(
    "-d",
    "--dir",
    required=False,
    type=click.Path(),
    show_default=True,
    default="",
    help="Path to the output dir.",
)
def pancore(modelpath, pcpath, based_on, dir):
    """Compare a model to a pan-core model."""
    model = rg.utility.io.load_model(modelpath, "cobra")
    pcmodel = rg.utility.io.load_model(pcpath, "cobra")
    report = rg.analysis.core_pan.compare_to_core_pan(model, pcmodel, based_on)
    report.save(dir)


@analyse.command()
@click.argument("modelpaths", nargs=-1, type=click.Path(exists=True))
@click.option(
    "--type",
    "-t",
    required=False,
    type=click.Choice(["sboterm"]),
    multiple=True,
    show_default=True,
    default=[],
    help="Type of comparison to be performed.",
)
@click.option(
    "--all",
    required=False,
    is_flag=True,
    show_default=True,
    default=False,
    help="Shortcut to run all comparisons. Overwrites input of type.",
)
@click.option(
    "-d",
    "--dir",
    required=False,
    type=click.Path(file_okay=False),
    show_default=True,
    default="comparison",
    help="Path to the output dir.",
)
@click.option(
    "-c",
    "--colour",
    required=False,
    type=str,
    default="Paired",
    show_default=True,
    help="Name of Matplotlib colour palette for the plot.",
)
def compare(modelpaths, type, all, dir, colour):
    """Compare models."""
    Path(dir).mkdir(parents=True, exist_ok=True)
    # compare SBOterms
    if all or "sboterm" in type:
        models = rg.utility.io.load_model(list(modelpaths), "libsbml")
        fig = rg.analysis.comparison.plot_rea_sbo_multiple(models, color_palette=colour)
        fig.savefig(Path(dir, "sboterm"), dpi=400, bbox_inches="tight")


# analyse growth
# --------------
@analyse.group()
def growth():
    """Analyse the growth under different conditions."""


# -> does it need to run without a media config as well?
@growth.command()
@click.argument("modelpaths", nargs=-1, type=click.Path(exists=True))
@click.option(
    "-m",
    "--media",
    required=True,
    type=click.Path(exists=True),
    help="Path to a media config file.",
)
@click.option(
    "-n",
    "--namespace",
    required=False,
    type=click.Choice(["BiGG"]),
    show_default=True,
    default="BiGG",
    help="Namespace to use for the model.",
)
@click.option(
    "-d",
    "--dir",
    required=False,
    type=click.Path(),
    show_default=True,
    default="",
    help="Path to the output dir.",
)
@click.option(
    "-c",
    "--colors",
    required=False,
    type=str,
    show_default=True,
    default="YlGn",
    help="Abbreviation of a matplotlib colour palette.",
)
def simulate(modelpaths, media, namespace, dir, colors):
    """Simulate the growth of the given model vs. media."""
    modelpaths = list(modelpaths)
    report = rg.analysis.growth.growth_analysis(
        modelpaths, media, namespace, retrieve="report"
    )
    report.save(to=dir, color_palette=colors)


@growth.command()
@click.argument("modelpath", type=click.Path(exists=True))
@click.option(
    "-m",
    "--media",
    required=True,
    type=click.Path(exists=True),
    help="Path to a media config file.",
)
@click.option(
    "-n",
    "--namespace",
    required=False,
    type=click.Choice(["BiGG"]),
    show_default=True,
    default="BiGG",
    help="Namespace to use for the model.",
)
@click.option(
    "-d",
    "--dir",
    required=False,
    type=click.Path(),
    show_default=True,
    default="",
    help="Path to the output dir.",
)
@click.option(
    "-c",
    "--colors",
    required=False,
    type=str,
    show_default=True,
    default="YlGn",
    help="Abbreviation of a matplotlib colour palette.",
)
def auxotrophies(modelpath, media, namespace, dir, colors):
    """Test for auxotrophies for the 20 proteinogenic amino acids."""
    medialist, supps = rg.classes.medium.load_media(media)
    model = rg.utility.io.load_model(modelpath, "cobra")
    report = rg.analysis.growth.test_auxotrophies(model, medialist, namespace)
    report.save(dir, colors)


@growth.command()
@click.argument("modelpath", type=click.Path(exists=True))
@click.option(
    "-e",
    "--element",
    type=str,
    required=True,
    help="Element to perform the source test with. Needs to be a valid chemical, elemental Symbol.",
)
@click.option(
    "-s",
    "--substances",
    type=str,
    required=False,
    show_default=True,
    default=None,
    multiple=True,
    help="Add substances from the database to test against. If none are given, tests against the whole database.",
)
@click.option(
    "-m",
    "--medium",
    required=False,
    type=str,
    show_default=True,
    default=None,
    help="Name of a medium from the database to use for the testing. If not given, uses the one from the model.",
)
@click.option(
    "-n",
    "--namespace",
    required=False,
    type=click.Choice(["BiGG"]),
    show_default=True,
    default="BiGG",
    help="Namespace to use for the model.",
)
@click.option(
    "-d",
    "--dir",
    required=False,
    type=click.Path(),
    default="",
    help="Path to the output dir.",
)
@click.option(
    "-c",
    "--colors",
    required=False,
    type=str,
    show_default=True,
    default="YlGn",
    help="Abbreviation of a matplotlib colour palette.",
)
def sources(modelpath, element, substances, medium, namespace, dir, colors):
    """Test growth on different sources for an element."""
    model = rg.utility.io.load_model(modelpath, "cobra")
    report = rg.analysis.growth.test_growth_with_source(
        model, element, substances, medium, namespace
    )
    report.save(dir, color_palette=colors)


@growth.command()
@click.argument("modelpath", type=click.Path(exists=True))
@click.option(
    "-o",
    "--objective",
    required=False,
    type=click.Choice(["flux", "medium", "exchanges"]),
    show_default=True,
    default="flux",
    help="Set the type of minimal medium to be calculated: minimal fluxes, minimal compounds or based on the model exchange reactions.",
)
@click.option(
    "-r",
    "--growth-rate",
    required=False,
    type=float,
    show_default=True,
    default=0.5,
    help="Minimal rate to be reached by the minimal medium. The smaller the rate the more expensive the calculation.",
)
@click.option(
    "-d",
    "--dir",
    required=False,
    type=click.Path(),
    show_default=True,
    default="",
    help="Path to the output dir.",
)
def minimal_medium(modelpath, objective, growth_rate, dir):
    """Calculate the minimal medium of a model.

    Either set the fluxes minimal, the number of compounds based on current medium or the number of exchange reactions.
    """

    model = rg.utility.io.load_model(modelpath, "cobra")
    medium = rg.analysis.growth.model_minimal_medium(model, objective, growth_rate)
    medium.export_to_file(dir=dir)
