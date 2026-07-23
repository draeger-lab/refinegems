SBOannotator with refineGEMs
============================

refineGEMs provides a connection to the `SBOannotator <https://github.com/draeger-lab/SBOannotator>`__\ :footcite:p:`Leonidou2023_sboann`
to annotate a model with SBO terms.

.. warning:: 

    The back-end of this function (:py:func:`~refinegems.utility.connections.run_SBOannotator`) is under heavy construction, but it should still work for the user.

The functionality can be called as follows:

.. code:: python
    :linenos:
    
    # model is the model loaded with libSBML
    from refinegems.utility.connections import run_SBOannotator

    model_with_sbo = run_SBOannotator(model)

.. footbibliography::
