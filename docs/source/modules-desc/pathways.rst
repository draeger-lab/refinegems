Addition of KEGG Pathways
=========================

The KEGG database holds information on metabolic pathways. 

Add KEGG pathways from reaction:
--------------------------------
You can use this to add KEGG pathways with the libSBML Groups plugin.

The workflow of the script is as follows:

1. Extraction of the KEGG reaction IDs from the annotations of your reactions
2. Identification, in which KEGG pathways these reactions occur
3. Addition of all KEGG pathways for a reaction with the biological qualifier ``OCCURS_IN`` to the annotations
4. Addition of all KEGG pathways as groups with references to the contained reactions as ``groups:member``

The main function for adding KEGG pathway groups is
:py:func:`~refinegems.curation.pathways.set_kegg_pathways`:

.. code:: python
    :linenos:
    
    from refinegems.curation.pathways import set_kegg_pathways
    from refinegems.utility.io import load_model, write_model_to_file

    model = load_model("path/to/model.xml", "libsbml")
    non_kegg_reactions = set_kegg_pathways(model)
    write_model_to_file(model, "path/to/model_with_pathways.xml")
