Usage as python module
======================

See :doc:`examples <examples>` to learn how to use refineGEMs.
Note that at this time most of the modules only make sense when you use the respective main functions:

.. autofunction:: refinegems.charges
    :noindex:
.. autofunction:: refinegems.polish
    :noindex:
.. autofunction:: refinegems.genecomp
    :noindex:
.. autofunction:: refinegems.modelseed
    :noindex:
.. autofunction:: refinegems.kegg_pathways
    :noindex:
.. autofunction:: refinegems.sboann
    :noindex:

The modules ``io``, ``cvterms`` and ``investigate`` provide functions that can be used by themselves.

``io`` 
--------
provides a couple of helper functions to load models, databases and parse gfffiles

.. automodule:: refinegems.io
    :members:
    :noindex:

``investigate``
--------
provides a couple of functions to get parameters of your model

.. automodule:: refinegems.investigate
    :members:
    :noindex:

``cvterms``
-----------
provides functions to work with cvterms

.. automodule:: refinegems.cvterms
    :members:
    :noindex:
