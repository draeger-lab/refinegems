Usage as python module
======================

Note that at this time most of the modules only make sense when you use the respective main functions:

.. autofunction:: refinegems.correct_charges
    :noindex:
.. autofunction:: refinegems.polish_carveme
    :noindex:
.. autofunction:: refinegems.genecomp
    :noindex:
.. autofunction:: refinegems.modelseed
    :noindex:
.. autofunction:: refinegems.kegg_pathways
    :noindex:
.. autofunction:: refinegems.sbo_annotation
    :noindex:

The modules ``load`` and ``test`` provide functions that can be used by themselves.

``load`` 
--------
provides a couple of helper functions to avoid remembering how cobra and libsbml actually load

.. automodule:: refinegems.load
    :members:
    :noindex:

``test``
--------
provides a couple of functions to get parameters of your model

.. automodule:: refinegems.test
    :members:
    :noindex:

``cvterms``
-----------
provides functions to work with cvterms

.. automodule:: refinegems.cvterms
    :members:
    :noindex:
