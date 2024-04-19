Usage
======

Access points
-------------

via command line/terminal
^^^^^^^^^^^^^^^^^^^^^^^^^

The main functionalities of refineGEMs can be accessed directly from the command line after a successfull installation.

For more information, call the following in your command line or terminal:

.. code-block:: bash

  refinegems --help


inside Python
^^^^^^^^^^^^^

Alternatively, refineGEMs can be directly imported into Python or a python script for full access of all modules, classes and more.

.. code-block:: python 

  import refinegems as rg

Or, e.g.:

.. code-block:: python

  from refinegems.analysis.growth import growth_analysis
  

How to use specific functionalities
-----------------------------------

TODO
- notebooks
- polish, gapfill, configs and more

See :doc:`examples <modules/examples>` to learn how to use refineGEMs.

Some additional remarks
-----------------------

.. warning:: 
    Using ``lab_strain=True`` has the following two requirements:
        
      1. The model already contains GeneProduct identifiers containing valid NCBI Protein/RefSeq identifiers.
          If there is no available data for the modeled organism in any database these identifiers can be added with 
          the pipeline described in :ref:`Pipeline: From genome sequence to draft model` before draft model creation.  
      2. Input of a FASTA file containing header lines similar to:
          >lcl|CP035291.1_prot_QCY37216.1_1 [gene=dnaA] [locus_tag=EQ029_00005] [protein=chromosomal replication initiator protein DnaA] [protein_id=QCY37216.1] [location=1..1356] [gbkey=CDS]
          Of the description part in the header line only locus_tag, protein and protein_id are important for :py:mod:`~refinegems.curation.polish`.
    


The media database
------------------

TODO