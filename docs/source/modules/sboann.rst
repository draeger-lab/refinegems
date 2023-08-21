SBOannotator with refineGEMs
============================

RefineGEMs offers access to the functionalities of `SBOannotator <https://github.com/draeger-lab/SBOannotator>`__\ :footcite:p:`Leonidou2023_sboann`. 

The ``sboann`` module is splitted into a lot of small functions which are all annotated, however when using it for SBO-Term annotation it only makes sense to run the function ``sbo_annotation``: 

.. autofunction:: refinegems.sboann.sbo_annotation
    :noindex:

.. code:: python
    :linenos:
    
    import refinegems as rg 
    model_sboann = rg.sboann.sbo_annotation(<path to your model>)
    rg.io.write_to_file(model_sboann, <path to modified model>)

If you use it from the refineGEMs toolbox with the config you can get a visualization of the SBO-Term distribution before and after the SBO-Term update.

.. footbibliography::
