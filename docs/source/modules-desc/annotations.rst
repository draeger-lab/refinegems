Annotation polishing
====================

Annotation polishing focuses on making existing model annotations more consistent and
MIRIAM-compliant. It is useful when a model already contains identifiers, but not if 
the URI patterns, qualifiers or annotation formatting are inconsistent.

The command-line entry point is:

.. code-block:: bash

    refinegems curate annotations MODEL.xml --outdir polished_annotations

Internally, this uses :py:func:`~refinegems.curation.miriam.polish_annotations`.
The function checks the annotation resources attached to model entities, generates
MIRIAM-compliant URI patterns and writes a table of invalid annotations, that could not be fixed, to the output
directory.

The option ``--new-pattern`` controls whether the newer CURIE pattern
``database_prefix:local_identifier`` is used when creating identifiers:

.. code-block:: bash

    refinegems curate annotations MODEL.xml --new-pattern True --outdir polished_annotations

When to use
-----------

Use annotation polishing before downstream workflows that depend on machine-readable
identifiers, for example pathway addition, gap-filling or external quality checks.
It can also be useful before comparing models, because normalised identifiers reduce
the number of differences caused only by formatting.
