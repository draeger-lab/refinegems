Correcting the biomass objective function
=========================================

To get feasible results from a model simulation, the biomass objective function (BOF) 
should add up to 1 mmol/gDW/h. However, this is often not the case, e.g. the BOF  of 
CarveMe does not add up to 1 mmol/gDW/h directly after the draft reconstruction.

Therefore, this needs to be corrected. 

The module :py:mod:`~refinegems.curation.biomass` provides the function :py:func:`~refinegems.curation.biomass.check_normalise_biomass` 
to check if the BOF of a model is consistent 
and if necessary adjusts/normalises the BOF to add up to 1 mmol/gDW/h. 

.. footbibliography::
