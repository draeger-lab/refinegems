Energy generating cycles (EGCs)
===============================

Energy generating cycles (EGCs) are cycles inside the model that are able to produce 
energy metabolites even when the model is cut of from any input. Since these cycles create 
energy basically out of nothing, they are biologically infeasible and are likely to distort 
the simulations, e.g. the growth simulation.

Therefore it is vital to identfy and resolve these cycles. The :py:mod:`~refinegems.classes.egcs` module 
provides classes to tackle these issues.

Identify EGCs
-------------

If the aim is to simply find the EGCs (e.g. for subsequent manual adjustment), 
The :py:class:`~refinegems.classes.egcs.EGCSolver` provides the function 
:py:meth:`~refinegems.classes.egcs.EGCSolver.find_egcs` to automatically find the egcs present in the model.
Its implementation is based on :footcite:t:`fritzemeier2017erroneous`:

    For each energy metabolite, the following steps are performed:

    1. A dissipation reaction for the energy metabolite is temporarily added to the model.
    2. The reaction bounds are limited to 1.0/-1.0 and all exchange reactions to the cytosol are closed.
    3. The dissipation reaction ist set as the model objective and the model is optimised.
    4. If the objective value is greater than the minimal growth threshold, at least one EGC for the corresponding energy metabolite is present.

This functionality can be accesed via

.. code-block:: bash 

    refinegems refine egcs 'MODELPATH' 

or 

.. code-block:: python

    from refinegems.classes.egcs import EGCSolver

    # model: cobra.Model
    solver = EGCSolver()
    results = solver.find_egcs(model)

By setting ``with_reacs`` to ``True``, the reactions that form the problematic cycle(s), meaning the one with non-zero flux values are additionally reported.


Solve EGCs
----------

Solving the EGCs is a challenging problem. Possible solution strategies include finding 
the least number of changes that need to be introduced into the model (purely mathematical approach)
or changing the reactions only according to biological evidence (purely biological approach).

Below is the list of currently available Solvers implemented in the toolbox. 

.. note::

    We are working on extending this list to allow for better, faster and/or more fitting 
    solutions for resolving EGCs.

Greedy approach
^^^^^^^^^^^^^^^

**class** :py:class:`~refinegems.classes.egcs.GreedyEGCSolver`

This approach is loosely based on :footcite:t:`fritzemeier2017erroneous`. 
It is a greedy, mathematical approach, that tries to resolve the EGCs by changing the reaction modes
of a single reaction each time and finding the smallest set of reactions that can resolve the biggest set of 
found EGCs. The modifications of the reactions include:

- Make reaction reversible
- Delete backward reaction
- Delete forward reaction
- Delete complete reaction

The algorithmn runs reasonably fast. However since it only considers single reaction modifications, 
some EGCs may remain unsolved. These will be reported at the end of the solving process.

This Solver can be accessed via

.. code-block:: bash 

    refinegems refine egcs -s greedy 'MODELPATH' 

or 

.. code-block:: python

    from refinegems.classes.egcs import GreedyEGCSolver

    # model: cobra.Model
    solver = GreedyEGCSolver()
    results = solver.solve_egcs(model)


.. footbibliography::