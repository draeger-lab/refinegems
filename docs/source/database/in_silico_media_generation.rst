From laboratory to *in silico* medium
==========================================

The `refineGEMs` toolbox allows user to add their own media definition to the workflows 
or to their models. Personal media definitions can be converted into a Medium object by calling 
the function :py:func:`~refinegems.classes.medium.read_external_medium` 
or added to the pipelines by using the `external_base` / `add_external` options in the media config.
The input for both is a path to a TSV file containing the medium information. In case of the functions, a 
command line options with prompt is also available. 

The format of the input TSV file and the conversion of a medium definition to an *in silico* medium
will be discussed in the following sections.

.. note::

   If you want to suggest to add your medium to the database, 
   contact the developers or comment in this `GitHub issue <https://github.com/draeger-lab/refinegems/issues/123>`_.



Medium file format of the toolbox
---------------------------------

The file format for a medium usable with the `refineGEMs` toolbox is a tab-separated file (TSV) 
in the following format (with example):

.. code::

   name  formula  flux  source   X  X  ....
   water H2O   10.0  water .....

This table is the substance table, containing all the substances and their information, that can be found in the medium.

The first non-comment line (more about that later) needs to be header. The TSV has to contain at least the first 
four columns as specified in :py:data:`~refinegems.classes.medium.REQUIRED_SUBSTANCE_ATTRIBUTES`. 
An arbitrary amount of additional column (represented by the *X* in the example above) 
can be added, however, only columns with the same as the ones listed in 
:py:data:`~refinegems.classes.medium.ALLOWED_DATABASE_LINKS` will be transferred into the Medium object. 
These are specifically the databases, that are covered (or will be soon) in the database and therefore kept for matching purposes.

In addition to the substances, more information about the medium can be added to the 
file as comment lines. To be correctly identified, these lines need to start with a `#`.

Additionally, one can add more information about the medium to the Medium object that will be created
by including information in the format:

.. code:: 

   # descriptor: value

If the descriptor is either *name*, *reference* or *description*, the information will 
be added as the *name*, *doi* and *description* attribute of the Medium object respectively.


How to get from a lab medium to the *in silico* one
---------------------------------------------------

To create your own *in silico* medium definition, follow the steps listed below:

1. Search papers containing medium definitions./ Search paper or provider information for a medium that could be 
   interesting for your organism.
2. If the paper contains already an *in silico* defintion: 

      i. Rearrange the data into the format described in :ref:`Medium file format of the toolbox`
      ii. Go to step 3.

   If not, create a table based on :ref:`Medium file format of the toolbox`:

      i. | Check which substances should be in the medium

         .. note::

            *Substance* means an indepandant chemical entity in the context of the medium, e.g. 
            the different ions of a salt are each one substance, even though the salt is one component 
            in the medium recipe, since chemically, the ions can react in reactions independant from one another, as 
            they dissolve in usually used aequeous solution.

         | If they already have an entry in the `refineGEMs` database, directly add the information to the TSV from there.
         | If they are not in the database, the information needs to be collected:
            
            a. The name of the substance needs to be unique, if it should later be added to the data, so choosing the 
               so the IUPAC name of the substance can be a good start. For better human readability, a trivial name can be added and the full name can be put in bracktes.
            b. Databases like ChEBI are a good source for the formula. Note, that in case of 
               multi-charge states, the charge state should correcpond to the name.
            c. The flux is dependant on the organism but is korrelated to the amount of the substance is the medium, 
               as the organism can at most take up as much substance as is in the medium 

            .. warning::

               explanation for conversion concentration -> flux missing
            
            d. The source describes where to substance originates from, e.g. if salt was 
               added to the medium, the substances Na+ and CL- need to be added to the medium, 
               but both originate from the same source, which is written down in the source column.
               Additionally, a substamce can originate from different source.
            e. Extent amd fill the database columns of xour choice, the more the better.

            .. hint::

               | ChEBI can be a good place to start looking for the substances and already provides some links to other databases.
               | Furthermore, MetaNetX is often well connected, so searching there first can speed up the process.

               When searching for the substance in the different databases, using different synomys e.g. found in ChEBI or the
               formula (with different charge / number of H-atoms) can improve the change of finding a hit.

            .. note::

               If you only want to use this medium for your model with known namespaces,
               you potentially could only at the information for that specific namespace.


      .. i. Collect all components (substances) in a table 
      .. ii. Add a column for the BiGG IDs of the components
      .. iii. | Search for all components in the BiGG database
      ..      | If component not in BiGG database:

      ..       a. Search component in PubChem & try to find component with similar name and formula in BiGG database
      ..       b. | Search component in KEGG database & provide KEGG identifier in parenthesis behind the component name
      ..          | (In this case the BiGG ID column remains empty.)

      .. iv. Add the Bigg ID name to the component column & Add brackets around the original component name
      .. v. | For salts create an extra table where each ion is listed separately
      ..    | -> For how the salts should be separated into their ions search the salt in PubChem and look at the 2D Structure.
      .. vi. Search each ion in the BiGG database & add if available (Otherwise just leave the column empty)
      .. vii. | Transfer the results for the salts into the original table:
         
      ..      |   Add the BiGG ID of each ion to the component containing the ion such that in each row only one ion is represented.
      ..      |   If duplicates occur only add the according ion to the first component containing it.

3. Check if the medium definition misses relevant components like water, iron, oxygen or carbondioxide.
   Especially trace components like iron are easily missed, as a miniscule amound is often enought to enable growth for bacteria.
   This amount can already be added by just using tap water instead of distilled water.
4. If relevant components are missing, search for reasonable explanations to add the according components or justify why 
   these are missing.
5. The medium is now ready to be used for growth simulation! :) 
