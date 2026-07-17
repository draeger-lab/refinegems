Available subsets
=================

The list of subsets available in the database can be found below. Each entry for a subset provides further information 
about the substances. The rules used to get the *in silico* subset definitions are the same as described for a medium 
(see Section :ref:`From laboratory to *in silico* medium`).

.. caution::

  The subset called ``DiReM`` is a pseudo subset as it holds no biological meaning when adding it to a medium and is 
  solely used in the :py:mod:`~refinegems.classes.egcs` module. It contains the metabolites needed for the dissipation 
  reactions, which in turn are used to check for the exsistence of energy generating cycles (EGCs).

.. toctree:: 
  :glob:
  
  subsets/*
