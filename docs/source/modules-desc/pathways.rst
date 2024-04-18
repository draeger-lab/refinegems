Addition of KEGG Pathways
=========================

The KEGG database holds information on metabolic pathways. You can use this module to add KEGG pathways with the libSBML Groups plugin.

The workflow of the script is as follows:
1. Extraction of the KEGG reaction IDs from the annotations of your reactions
2. Identification, in which KEGG pathways these reactions occur
3. Addition of all KEGG pathways for a reaction with the biological qualifier ``OCCURS_IN`` to the annotations
4. Addition of all KEGG pathways as groups with references to the contained reactions as ``groups:member``

The only function that you will need to access is ``kegg_pathways``:

.. autofunction:: refinegems.pathways.kegg_pathways
    :noindex:

.. code:: python
    :linenos:
    
    import refinegems as rg 
    model_pathway_groups, non_kegg_reactions = rg.pathways.kegg_pathways(<path to your model>)
    rg.io.write_to_file(model_pathway_groups, <path to modified model>)