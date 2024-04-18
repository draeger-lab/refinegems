Growth simulation with refineGEMs
=================================

Interpreting the results
------------------------
Outputs a table with the column headers:

- essential = metabolites not present in the medium but necessary for growth
- missing = all exchanges missing in the model but given in medium


Implementation
--------------

Growth rates and, thus, doubling times can be determined with Flux Balance Analysis (FBA). 
RefineGEMs uses a COBRApy based implementation that is able to add metabolites to media definitions until growth is obtained. 
The pseudocode is shown below.

.. sourcecode:: python
  :linenos:
  :caption: Pseudocode representation of the algorithm implemented for growth simulation.

  model  # a model loaded with COBRApy
  medium # a Medium object, loaded from ex- or intern 

  ex_medium = medium.medium_to_model() # export to model

  if supplementation == True:
    new_medium = ex_medium.find_additives_to_enable_growth(params)
  else:
    new_medium = ex_medium

  model.medium = new_medium
  report = growth_simulation(model)

  return report

Options for supplementation currently include ``std`` for using the standard medium uptake or ``min`` for using the minimum medium uptake for supplementation.
Both are calculated by the ``get_uptake`` function of the module.

Usable Media
------------
The media to be used with the ``growth`` module via a ``config.yaml`` file.
More information about which media are provided by the database and how to generate *in silico* media can be found in Section :ref:`Media & Subsets`.


