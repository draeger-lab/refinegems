Growth simulation with refineGEMs
=================================

Interpreting the results
------------------------
Outputs a table with the column headers:

- essential = metabolites not present in the medium but necessary for growth
- missing = all exchanges missing in the model but given in medium


Implementation
--------------

Growth rates and, thus, doubling times can be determined with Flux Balance Analysis (FBA). RefineGEMs uses a COBRApy based implementation that adds metabolites one-by-one to custom media definitions until growth is obtained. The pseudocode is shown below.

.. image:: ../images/growth_algorithm.png
  :align: center
  :width: 400
  :alt: Pseudocode representation of the algorithm implemented for growth simulation.

There is a flag called basis which can be set to either ``default_uptake`` or ``minimal_uptake``. You can decide from which uptake you want to fill your medium of interest when looking for missing metabolites. Either the ``default_uptake`` which is the uptake that the model has when no specific medium is set or the ``minimal_uptake`` which is the uptake resulting from COBRApy's minimal_medium optimization.

Usable Media
------------
The media to be used with the ``growth`` module via a ``config.yaml`` file.
More information about which media are provided by the database and how to generate *in silico* media can be found in Section :ref:`Media`.


