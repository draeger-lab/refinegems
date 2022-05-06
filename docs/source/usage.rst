Usage as python module
======================

Note that at this time most of the modules only make sense when you use the respective main functions:

.. autofunction:: refinegems.charge
.. autofunction:: refinegems.polish_carveme
.. autofunction:: refinegems.genecomp
.. autofunction:: refinegems.modelseed
.. autofunction:: refinegems.kegg_pathways
.. autofunction:: refinegems.sbo_annotation

The modules ``load`` and ``test`` provide functions that can be used by themselves.

``load`` 
--------
provides a couple of helper functions to avoid remembering how cobra and libsbml actually load

.. automodule:: refinegems.load
    :members:

``test``
--------
provides a couple of functions to get parameters of your model

.. automodule:: refinegems.test
    :members:

