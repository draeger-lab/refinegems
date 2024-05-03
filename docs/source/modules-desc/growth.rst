Growth simulation with refineGEMs
=================================

The :py:mod:`~refinegems.analysis.growth` submodule provides a number of functionalities
to analyse and simulate the growth of a model including:

1. actual growth simulation on different media
2. auxothrophy simulation
3. simulating growth on different elemental sources, e.g. a C-source test
4. model the minimal medium for a model

Growth 
------

Using :py:func:`~refinegems.analysis.growth.growth_analysis`, the growth between any number of models and media can be simulated.
For easier use, the media can be describes using a ``config_media.yaml`` file, which in turn can be used as input for the function.
The media configuration file allows for many kinds of manipulation of the media from the in-build database, 
as well as adding additional media from external sources to the list of media to be tested.

.. hint::

  More information about the in-build media can be found under :ref:`Media & Subsets`.

Furthermore, a supplementation mode can be set to ensure the model grows on the given medium.
Supplemented compounds are reported at the end to ensure total transparency.

.. note::

  Currently implemented supplementation mode (in addition to none) include:

  - std/standard: supplementation based on standard, non-zero flux medium compounds
  - min/minimal: supplementation from the minimal medium the model is able to grow on

The function return the report, plot or both depending on the input for ``retrieve``.
The report contains the growth value and doubling times, supplemented compounds and 
compounds for which no exchange reactions was found in the model. The report also proves functions for 
further visualisation of these values,including a tabulary and graphical display.

Auxothrophy 
-----------

The function :py:func:`~refinegems.analysis.growth.test_auxotrophies` tests the auxothrophies
for the 20 proteinogenic amino acids for a model, a set of media and supplement modes.

Iteratively for each amino acids, a sink reactions for that particular amino acids is temporarily added
to the model and set as the objective value. In the case the model has an exchange reaction for the 
amino acids, it is disabled to ensure, the only possible way for the model to have the amino acids is
by producing it itself. The objective values for the added sink reaction are reported. 
Values higher than the growth threshold signify the ability of the model to produce the corresponding
amino acids by itself (the model is auxotroph for the amino acid).

The function returns a report, which in turn can be e.g. visualised as a heatmap.

Below, an examplary heatmap based on the Master Thesis of Carolin Brune (...) for 
*Klebsiella pneumoniae* MD01 can be seen. The model is able to produce all amino acids 
and therefore is autotroph for all amino acids.

.. image:: ../images/MT_CB_Kpneu_auxo.png


Source test
-----------

Using the :py:func:`~refinegems.analysis.growth.test_growth_with_source` function,
the user can simulate the growth of a model on a medium while varying the source of a given
element, e.g. carbon. Using the parameter ``substances``, alternative sources can be set.
If none are given, the program tests for all available sources in the in-build database for the given element.

Iteratively, the function changes all metabolites containing the given element with the next in list.
The model is optimised and the objective values for the BOF ist returned for all the different compounds.

.. note::

  For sensible test results, the starting medium should provide all other element necessary
  for the growth of the model as other compounds. Only under these conditions can the results of the 
  source test savely discussed in the context of changung a single source only.

Below you can see the exempkary results for a *Klebsiella pneumoniae* model with a supplemented M9 medium.
The supplemented M9 medium contains all trace elements the bacterium needs to grow e.g. iron and a single 
nitrogen source that is swapped out during the test. Since the original nitrogen source was ammonia, swapping 
it out does not eliminate another element from the medium.

.. image:: ../images/N-source.png 


Minimal medium
--------------

Using the :py:func:`~refinegems.analysis.growth.model_minimal_medium`, a minimal medium (composition)
can be calculated.

Since the description of a minimal medium can vary, refineGEMs provides different 
ways to calculate one, set using the parameter ``objective``:

- flux: find the minimal fluxes for the current medium 
- medium: find the minimal number of compounds based on the current medium
- exchanges: find the minimal number of compounds for a medium based on the exchange reactions in the model

.. note:: 

  The function always returns a single solution, but there may be more than one solution, especially
  for the third case.


