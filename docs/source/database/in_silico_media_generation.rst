From laboratory to *in silico* medium
==========================================

The `refineGEMs` toolbox allows users to add their own media definition to the workflows 
or to their models. Personal media definitions can be converted into a Medium object by calling 
the function :py:func:`~refinegems.classes.medium.load_external_medium`
or added to the pipelines by using the ``external_base`` / ``add_external`` options in the media config.
The input for both is a path to a tab-separated (TSV) file containing the medium information.

.. 
   @ASK How to reformulate pipelines here? That is not really fitting as refineGEMs does not contain pipelines.

The format of the input TSV file and the conversion of a medium definition to an *in silico* medium
will be discussed in the following sections.

.. hint::

   Please contact the developers of leave a comment in this 
   `GitHub issue <https://github.com/draeger-lab/refinegems/issues/123>`_ for any suggestions on new media for the 
   internal database.


Medium file format of the toolbox
---------------------------------

The file format for a medium usable with the `refineGEMs` toolbox is a TSV
in the following format (with example):

.. literalinclude:: ../../../src/refinegems/example/example_inputs/custom_medium_substance_table.tsv

This table is the substance table, containing all the substances, that can be found in the medium, and their corresponding information.

The first non-comment line needs to be the header. The TSV has to contain at least the first 
four columns as specified in :py:data:`~refinegems.classes.medium.REQUIRED_SUBSTANCE_ATTRIBUTES`. 
An arbitrary amount of additional columns can be added. However, only columns with a name listed in 
:py:data:`~refinegems.classes.medium.ALLOWED_DATABASE_LINKS` will be transferred into the Medium object. 
These are specifically the databases, that are covered (or soon to be covered) in the database and are therefore kept for mapping purposes.

In addition to the substances, more information about the medium can be added to the 
file as comment lines. To be correctly identified, these lines need to start with a `#`.
If the descriptor is either *name*, *reference* or *description*, the information will 
be added as the *name*, *doi* and *description* attribute of the Medium object, respectively.


How to get from a lab medium to the *in silico* one
---------------------------------------------------

To create your own *in silico* medium definition, follow the steps listed below:

1. Search papers containing medium definitions./ Search paper or provider information for a medium that could be 
   interesting for your organism.
2. If the paper already contains an *in silico* definition: 

      i. Rearrange the data into the format described in :ref:`Medium file format of the toolbox`
      ii. Go to step 3.

   If not, create a table based on :ref:`Medium file format of the toolbox`:

      i. | Check which substances should be in the medium.

         .. note::

            *Substance* means an independant chemical entity in the context of the medium, e.g., 
            an ion in a salt. Each type of ion in a salt is one substance, even though the salt is one component 
            in the medium recipe, since chemically, they dissolve in aqueous solution. Thus, they can 
            react independently from one another.

         | If the substances already have an entry in ``refineGEMs``' database, directly add the information to the TSV.
         | If they are not in the database, the information needs to be collected:
            
            a. The name of the substance needs to be unique if it should later be added to the database. Hence, 
               choosing the IUPAC name of the substance can be a good start. For better human readability, a 
               trivial name can be added, and the full name can be put in brackets instead.
            b. Databases like ChEBI are good sources for the formula. Note that in the case of 
               multi-charge states, the charge state should correspond to the name and formula.
            c. The flux is dependent on the organism but is correlated to the amount of the substance in the medium, 
               as the organism can at most take up as much substance as is in the medium.

            .. caution::

               Converting concentrations to fluxes is difficult.
               Getting the right flux values for your organism directly from a laboratory would be the best option.
            
            d. The source describes where the substance originates from, e.g. if salt was 
               added to the medium, the substances Na+ and CL- need to be added to the medium, 
               but both originate from the same source, which is written down in the source column.
               Additionally, a substance can originate from different sources.
            e. Extend and fill the database columns of your choice; the more, the better.

            .. hint::

               | ChEBI can be a good place to start looking for the substances and often provides links to other databases.
               | Furthermore, MetaNetX is often well-connected, so searching there first can speed up the process.

               When searching for the substance in the different databases, using different synonyms, e.g. found in ChEBI or the
               formula (with different charge / number of H-atoms), can improve the chance of finding a hit.

            .. note::

               If you want to use this medium only for your model with known namespaces, 
               you could potentially add the information only for that specific namespace.

3. Check if the medium definition lacks relevant components like water, iron, oxygen or carbon dioxide.
   Trace components like iron are easily missed, as a miniscule amount is often enough to enable growth of bacteria.
   This amount can already be added by just using tap water instead of distilled water.
4. If relevant components are missing, search for reasonable explanations to add the corresponding components or justify why 
   these are missing.
5. The medium is now ready to be used for growth simulation! :)
