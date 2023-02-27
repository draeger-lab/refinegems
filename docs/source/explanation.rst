Ideas behind some implementations
=================================

Below you will find some arguing behind the implementation choices we made.

Gap filling
-----------

The `gapfill` module can be used either with KEGG were you only need the KEGG organism ID or with BioCyc. For the gap filling using BioCyc you will need to create an  `account <https://biocyc.org/>`. Then you will search for the strain of your organism. Then you click on Tools and select Special Smart Tables. There you can download the tables "all genes of <organism>" and "all reactions of <organism>". Remove all columns except Gene Name and then choose a transform (reactions of gene), then add property (Accession-2) and after that delete the Gene Name files. The choose Accession-2, and use the filter function to delete all empty rows. The Export to Spreadsheet File then choose frame ids (for gene table).

For reaction: 'Reaction' (present) 'Reactants of reaction' (transform) 'Products of reaction' (transform) 'EC-Number' (present) 'Reaction-Direction' (property) 'Spontaneous?' (property)

For metabolite: Use MetaCyc database -> all compounds
'Compound' (present) 'Object ID' (property) 'Chemical Formula' (property) 'InChI-Key' (property) 'ChEBI' (property - database links - ChEBI)
Export with common names

Running times: Around 2h for KEGG
Around 45 mins for BioCyc 

Growth simulation
-----------------

Growth rates and thus doubling times can be determined with Flux Balance Analysis (FBA). RefineGEMs uses a COBRApy based implementation that adds metabolites one-by-one to custom media definitions until growth is obtained. The pseudocode is shown below.

.. image:: images/growth_algorithm.png
  :width: 400
  :alt: Pseudocode representation of the algorithm implemented for growth simulation.

There is a flag called basis which can be set to either default_uptake or minimal_uptake. You can decide from which uptake you want to fill your medium of interest when looking for missing metabolites. Either the default_uptake which is the uptake that the model has when no specific medium is set or the minimal_uptake which is the uptake resulting from cobrapys minimal_medium optimization.