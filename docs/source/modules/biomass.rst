Correcting the biomass objective function
=========================================

The biomass objective function (BOF) of CarveMe does not add up to 1 mmol/gDW/h directly after the draft reconstruction.
However, to get feasible results from a model this reaction should add up to 1 mmol/gDW/h.
Thus, the module `biomass` was created to check if the BOF of a model is consistent and if necessary the BOF is adjusted 
to add up to 1 mmol/gDW/h.

.. footbibliography::
