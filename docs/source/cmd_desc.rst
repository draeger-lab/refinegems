Command Line access points
==========================

After a successful installation, ``refineGEMs`` can be accessed via the command line
from inside the Python environment it was installed in using:

.. code-block:: bash

    refinegems [OPTIONS] COMMAND [ARGS]

The following command groups are available:

- ``analyse``: Provides commands for analysing one or more models.
- ``curate``: Provides commands for polishing model annotations and model structure.
- ``database``: Provides commands for handling the internal database.
- ``media``: Provides commands for interacting with the media database.
- ``refine``: Provides commands for refining a model, including biomass correction,
  charge correction, gap-filling, reaction direction correction, EGC checks and annotation.
- ``setup``: Provides commands for downloading files and setting up supporting data.

General Options
---------------

- ``--help``/``-h``: Call the help page of the command.
- ``--version``: Show the version and exit.

----

refinegems analyse
------------------

Analyse a model by testing growth, finding EGCs, comparing models and reporting
overall model statistics.

.. code:: bash

  refinegems analyse compare [MODELPATHS]

Compare multiple models.

Options:

- ``--type/-t <sboterm>``: Type of comparison to perform. Multiple can be added.
- ``--all``: Perform all comparisons. Overwrites previous option.
- ``--dir/-d``: Path to the output directory.
- ``--colour/-c``: Name of a Matplotlib colour palette for the plot.

.. code:: bash

  refinegems analyse memote [MODELPATH]

Perform a memote analysis on a model.

Options:

- ``--score-only/-s``: Return only the final score of the analysis.
- ``--file/-f``: Name/path to save output to. Only relevant if ``-s`` was not set.

.. code:: bash

  refinegems analyse pancore [MODELPATH] [PCPATH]

Compare the pan-core information content of a model to a pan-core model.

Options:

- ``--based-on/-b``: How to compare the models. Currently supported: ``id``.
- ``--dir/-d``: Path to the output directory.

.. code:: bash

  refinegems analyse pathways [MODELPATH]

Analyse the KEGG pathways contained in a model.

Options:

- ``--dir/-d``: Path to the output directory.
- ``--colors/-c``: Abbreviation of a Matplotlib colour palette for the graphics.

.. code:: bash

  refinegems analyse stats [MODELPATH]

Perform a statistical analysis of the model.

Options:

- ``--dir/-d``: Path to the output directory.
- ``--colors/-c``: Abbreviation of a Matplotlib colour palette for the graphics.

refinegems analyse growth
^^^^^^^^^^^^^^^^^^^^^^^^^

Analyse the growth under different conditions.

.. code:: bash

  refinegems analyse growth auxotrophies [MODELPATH]

For a given set of media, simulate the amino acid auxotrophies of the model.

Options:

- ``--media/-m``: Path to a media config file. Required.
- ``--namespace/-n``: Namespace of the model.
- ``--colors/-c``: Abbreviation of a Matplotlib colour palette for the graphics.
- ``--dir/-d``: Path to the output directory.

.. code:: bash

  refinegems analyse growth minimal-medium [MODELPATH]

Calculate the minimal medium of a model. It can either be calculated by minimising
the fluxes of the current medium (``flux``), finding the minimal number of compounds
needed for growth based on the current medium (``medium``) or finding the minimal
number of compounds based on the available exchange reactions (``exchanges``).

Options:

- ``--objective/-o``: One of ``flux``, ``medium`` or ``exchanges``.
- ``--growth-rate/-r``: Minimal growth rate that should be reached on the minimal medium.
- ``--dir/-d``: Path to the output directory.

.. code:: bash

  refinegems analyse growth simulate [MODELPATHS]

Simulate the growth of one or multiple models on one or more media.

Options:

- ``--media/-m``: Path to a media config file. Required.
- ``--namespace/-n``: Namespace of the model.
- ``--colors/-c``: Abbreviation of a Matplotlib colour palette for the graphics.
- ``--dir/-d``: Path to the output directory.

.. code:: bash

  refinegems analyse growth sources [MODELPATH]

Simulate growth for different sources. When choosing the starting medium,
it is useful to provide at least one additional source of every other necessary
element separately from the one that is tested.

Options:

- ``--element/-e``: Element to perform the source test for. Must be a valid chemical symbol.
- ``--substances/-s``: Substances for substitution. Multiple can be given. If none are given, all options found in the database will be used.
- ``--medium/-m``: Medium abbreviation of a medium in the database, if the medium currently saved with the model should be substituted before testing.
- ``--namespace/-n``: Namespace of the model.
- ``--colors/-c``: Abbreviation of a Matplotlib colour palette for the graphics.
- ``--dir/-d``: Path to the output directory.

refinegems curate
-----------------

Curate or polish model annotations and model structure.

.. code:: bash

  refinegems curate annotations [MODELPATH]

Polish annotations by changing qualifiers and annotations to be MIRIAM-compliant.

Options:

- ``--new-pattern/-n``: Use the ``database_prefix:local_identifier`` CURIE pattern.
- ``--outdir/-o``: Path to a directory to write the output files to.

.. code:: bash

  refinegems curate model [MODELPATH]

Run the model-polishing workflow for a model created by an automatic reconstruction
pipeline. This extends annotations, cleans note fields and changes qualifiers and
annotations to be MIRIAM-compliant.

Options:

- ``--id-db/-i``: Main database where identifiers in the model come from.
- ``--mtf/--mapping-tbl-file/-m``: Path to a gene product mapping table.
- ``--lab-strain/-l``: Set to ``True`` for strains with only homolog mappings.
- ``--kid/--kegg-organism-id/-k``: KEGG organism identifier, if available.
- ``--outdir/-o``: Path to a directory to write the output files to.
- ``--gff-paths/-g``: Path(s) to GFF file(s). Used when no mapping table is provided.
- ``--email/-e``: E-mail for NCBI queries. Used when no mapping table is provided.
- ``--lt/--contains-locus-tags/-t``: Set to ``True`` if the model has locus tags within the label tag.

refinegems database
-------------------

Access and curate the internal database.

.. code:: bash

  refinegems database initialise

Initialise or update the internal database.

.. code:: bash

  refinegems database add-namespace [DATABASENAME]

Add or update tables for additional namespaces/databases in the internal database.

Options:

- ``--chunksize/-c``: Size in kB of data to download per chunk, if a download is required.

.. code:: bash

  refinegems database reset

Reset the database by removing additionally downloaded tables.

refinegems media
----------------

Access the media part of the database.

.. code:: bash

  refinegems media info

Retrieve information about the media database.

Options:

- ``--list``: List the available media.

.. code:: bash

  refinegems media export

Export media from the internal media database based on a medium name, a list of
medium names, ``all`` or a media configuration file.

Options:

- ``--media-names-or-config/--mn-cfg/-m``: Name(s) of media to export or path to a media configuration file.
- ``--type/-t``: Export file type. One of ``tsv``, ``csv``, ``docs`` or ``rst``.
- ``--flavour/-f``: Export flavour for ``tsv`` and ``csv``. One of ``substance_table`` or ``carveme_mimic``.
- ``--single-file/-s``: Export all media in one file. Only viable for ``carveme_mimic``.
- ``--no-flux/-n``: Remove the flux column from ``tsv``/``csv`` exports.
- ``--dir/-d``: Output directory.
- ``--max-widths/-w``: Maximal table width for documentation tables.

refinegems refine
-----------------

Refine a model. This includes steps like biomass correction, charge correction,
gap-filling, reaction direction correction, EGC checks, SBO annotation and pathway
annotation.

.. code:: bash

  refinegems refine automated-gapfill [ALG] [MODELPATH]

Find gaps in a model and optionally try to fill them. ``ALG`` must be one of
``KEGG``, ``BioCyc`` or ``Gene``.

General options:

- ``--outdir/-o``: Path to a directory to write the output to.
- ``--fill/-f``: Try to fill the gaps in the model.
- ``--formula-check/-fc``: One of ``none``, ``existence``, ``wildcard`` or ``strict``.
- ``--no-dna``: Exclude DNA reactions (name-based) from being added to the model.
- ``--no-rna``: Exclude RNA reactions (name-based) from being added to the model.
- ``--idprefix/-p``: Prefix for random IDs if an ID does not exist for the given namespace.
- ``--namespace/-n``: Namespace used in the model.
- ``--report/-r``: Save statistics and genes/reactions for manual curation.

KEGG required options:

- ``--orgid``: KEGG organism ID.

BioCyc required options:

- ``--genetable/-gt``: Path to the BioCyc gene SmartTable.
- ``--reactable/-rt``: Path to the BioCyc reaction SmartTable.
- ``--gff-bc``: Path to the GFF.

Gene required options:

- ``--gff-g``: Path to the GFF.

Gene optional options:

- ``--prot-prefix``: Prefix for pseudo-protein IDs.
- ``--mail``: Mail address for NCBI requests.
- ``--check-ncbi/-ncbi``: Enable searching protein IDs in NCBI.
- ``--db-type``: Database to search against. One of ``swissprot`` or ``user``.
- ``--fasta``: Path to the protein FASTA of the model.
- ``--dmnd-db``: Path to the DIAMOND database.
- ``--db-mapping/-db-map``: Path to the SwissProt or user-defined mapping file.
- ``--sensitivity/-s``: DIAMOND sensitivity mode.
- ``--cov``: Coverage value passed to DIAMOND.
- ``--pid``: Percentage identity threshold for filtering DIAMOND results.
- ``--threshold-add-reacs``: Maximal number of reactions for one EC number mapping for it to be added.
- ``--threads/-t``: Number of DIAMOND threads.

Constraints:

- ``--mail`` is required if ``--check-ncbi`` is set.
- If one of ``--db-type``, ``--fasta``, ``--dmnd-db`` or ``--db-mapping`` is set, all four need to be set.

.. code:: bash

  refinegems refine biomass [MODELPATH]

Normalise the biomass objective function(s) of a model to improve model consistency.

Options:

- ``--cycles/-c``: Maximal number of normalisation cycles.
- ``--outfile/-o``: Filename to save the updated model under.

.. code:: bash

  refinegems refine charges [MODELPATH]

Compare the charges in a model to the ModelSEED database and adjust them accordingly,
if necessary.

Options:

- ``--dir/-d``: Path to the output directory.

.. code:: bash

  refinegems refine direction [MODELPATH] [DATA]

Check and, if necessary, correct the direction of reactions in a model.

Options:

- ``--dir/-d``: Path to the output directory.

.. code:: bash

  refinegems refine egcs [MODELPATH]

Identify problematic energy generating cycles (EGCs) in the model and optionally try
to resolve them.

Options:

- ``--solver/-s``: Solver to use. Currently available: ``greedy``.
- ``--namespace/-n``: Namespace of the model.
- ``--compartment/-c``: Compartments to check, separated by comma. Defaults to ``c,e``.
- ``--outfile/-o``: File to save the updated model to, if a solver has been set.

refinegems refine annot
^^^^^^^^^^^^^^^^^^^^^^^

Add annotations to your model.

.. code:: bash

  refinegems refine annot sboterms [MODELPATH]

Call SBOannotator on a model to enhance/add SBO terms to the annotations.

Options:

- ``--dir/-d``: Path to the output directory.

.. code:: bash

  refinegems refine annot pathways [MODELPATH]

Add KEGG pathways as group entities to the model.

Options:

- ``--dir/-d``: Path to the output directory.

refinegems setup
----------------

Set up tools, folder structure and data for running the program.

.. code:: bash

  refinegems setup config

Download a configuration file for a specific functionality of the toolbox.

Options:

- ``--filename/-f``: Name or path of a file to save the config under.
- ``--type/-t``: Type of config to download. One of ``media`` or ``refinegems``.

.. code:: bash

  refinegems setup data [DOWNLOADTYPE]

Download files needed for a given functionality of the toolbox.

Current options include:

- ``SwissProt_gapfill``: Download the SwissProt files for gap-filling.

Options:

- ``--dir/-d``: Path to directory to save the downloaded files to.
- ``--chunksize/-c``: Size of the chunk to download in kB.

.. code:: bash

  refinegems setup geneproduct-mapping-table [MODELPATH]

Generate an ID mapping file for gene products.

Options:

- ``--gff-paths/-g``: Path(s) to GFF file(s). Allowed GFF formats are RefSeq, NCBI and Prokka.
- ``--email/-e``: E-mail for NCBI queries.
- ``--lt/--contains-locus-tags/-t``: Set to ``True`` if the model has locus tags within the label tag.
- ``--outdir/-o``: Path to a directory to write the output files to.

.. code:: bash

  refinegems setup build-pancore [MODELS]

Using the given models, construct a pan-core model.

Options:

- ``--based-on/-o``: Option on how to compare the models. Currently supported: ``id``.
- ``--name/-n``: Set the name of the constructed pan-core model.
- ``--keep-genes/-g``: Keep genes in the pan-core model.
- ``--rcomp/--resolve-compartments``: Try to standardise the model's compartment names.
- ``--dir/-d``: Path to the output directory.
