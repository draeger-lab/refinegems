SBOannotator with refineGEMs
============================

refineGEMs provides a connection to the `SBOannotator <https://github.com/draeger-lab/SBOannotator>`__\ :footcite:p:`Leonidou2023_sboann`
to annotate a model with SBO terms.

The functionality can be called as follows via Python:

.. code:: python
    :linenos:
    
    # model is the model loaded with libSBML
    from refinegems.utility.connections import run_SBOannotator

    model_with_sbo = run_SBOannotator(model)

or via the command line:

.. code:: bash

    refinegems refine annot sboterms MODEL.xml --dir out_dir

.. footbibliography::
