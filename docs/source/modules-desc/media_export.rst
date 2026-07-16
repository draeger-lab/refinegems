Media export
============

The internal media database can be exported into files for inspection, reuse or for
generating documentation tables.

Export all media
----------------

To export all media definitions as TSV files:

.. code-block:: bash

    refinegems media export --media-names-or-config all --type tsv --dir media_export

The options can also be shortened to:

.. code-block:: bash

    refinegems media export -m all -t tsv -d media_export

Export selected media
---------------------

Multiple media names can be passed by repeating ``-m``:

.. code-block:: bash

    refinegems media export -m LB -m M9 --type csv --dir selected_media

Export from a media configuration file
--------------------------------------

The same command can use a media configuration file instead of explicit medium names:

.. code-block:: bash

    refinegems media export -m media_config.yml --type tsv --dir configured_media

Export formats
--------------

The ``--type`` option controls the output format:

- ``tsv`` or ``csv`` for tabular files
- ``docs`` or ``rst`` for documentation-oriented tables

For ``tsv`` and ``csv`` exports, ``--flavour`` controls the table layout. The default
``substance_table`` export keeps the refineGEMs substance table format. The
``carveme_mimic`` export creates files closer to CarveMe-style media definitions.

Example output
--------------

For example, exporting M9 as ``substance_table`` starts with:

.. code-block:: text

    name	formula	flux	source	db_id	db_type
    Ammonia	H3N		NH4Cl	C00014	KEGG
    Ammonia	H3N		NH4Cl	C01342	KEGG
    Ammonia	H3N		NH4Cl	MNXM729302	MetaNetX
    Ammonia	H3N		NH4Cl	nh3	BiGG

The same medium exported as ``carveme_mimic`` starts with:

.. code-block:: text

    medium	description	compound	name	flux
    M9	M9 minimal medium	nh3	Ammonia	10.0
    M9	M9 minimal medium	nh4	Ammonia	10.0
    M9	M9 minimal medium	ca2	Calcium(II) cation [Ca(2+)]	10.0
    M9	M9 minimal medium	cl	Chloride [Cl(-)]	10.0
