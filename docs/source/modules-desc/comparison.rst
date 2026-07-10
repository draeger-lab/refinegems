Pan-core analysis and model comparison
======================================

``refineGEMs`` contains functions for working with pan-core models and for comparing model content. 
These functions are useful when curating a set of related strain models or when checking how one model 
differs from a reference or pan-core model.

.. note::

    The pan-core module is not fully developed yet. The current implementation
    should be treated as a starting point and may be extended or revised in a
    future update.

Pan-core model construction
---------------------------

A pan-core model can be constructed from multiple input models with:

.. code-block:: bash

    refinegems setup build-pancore MODEL_A.xml MODEL_B.xml MODEL_C.xml --name pan-core-model --dir pancore

By default, genes are removed from the generated pan-core model. Use ``--keep-genes``
if the gene layer should be retained:

.. code-block:: bash

    refinegems setup build-pancore MODEL_A.xml MODEL_B.xml --keep-genes --dir pancore

The Python API for this workflow is
:py:func:`~refinegems.analysis.core_pan.generate_core_pan_model`.

Comparing a model to a pan-core model
-------------------------------------

Once a pan-core model exists, another model can be compared against it:

.. code-block:: bash

    refinegems analyse pancore MODEL.xml PANCORE.xml --dir pancore_comparison

The corresponding Python function is
:py:func:`~refinegems.analysis.core_pan.compare_to_core_pan`.

Comparing SBO annotations
-------------------------

The command-line comparison currently supports comparing SBO term distributions across
multiple models:

.. code-block:: bash

    refinegems analyse compare MODEL_A.xml MODEL_B.xml --type sboterm --dir comparison

To run all available comparison plots, use:

.. code-block:: bash

    refinegems analyse compare MODEL_A.xml MODEL_B.xml --all --dir comparison

The Python API behind this workflow is available in
:py:mod:`~refinegems.analysis.comparison`, especially
:py:func:`~refinegems.analysis.comparison.sbo_terms`.
