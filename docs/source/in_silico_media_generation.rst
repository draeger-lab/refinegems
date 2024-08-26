From laboratory to *in silico* medium
======================================

.. warning:: 
    | *Will be deprecated from version 2.0.0 onwards.*
    | Using the new setup, the generation of *in silico* media changed slightly. 
    | Please refer to the newer documentation page after the version change or when using `SPECIMEN <https://github.com/draeger-lab/SPECIMEN>`__.


.. hint:: 
   If you want to use the medium with ``refineGEMs.growth`` add the definition to the database schema ``sbo_media_db.sql`` 
   in the folder *refinegems/database* in the downloaded repository. To update the database with the newly added table just 
   delete the file ``data.db`` in the same folder and run refineGEMs.

1. Search papers containing medium definitions./ Search paper or provider information for a medium that could be 
   interesting for your organism.
2. | If the paper contains already an *in silico* defintion: Go to step 3.
   | If not:

      i. Collect all components in a table 
      ii. Add a column for the BiGG IDs of the components
      iii. | Search for all components in the BiGG database
           | If component not in BiGG database:

            a. Search component in PubChem & try to find component with similar name and formula in BiGG database
            b. | Search component in KEGG database & provide KEGG identifier in parenthesis behind the component name
               | (In this case the BiGG ID column remains empty.)

      iv. Add the Bigg ID name to the component column & Add brackets around the original component name
      v. | For salts create an extra table where each ion is listed separately
         | -> For how the salts should be separated into their ions search the salt in PubChem and look at the 2D Structure.
      vi. Search each ion in the BiGG database & add if available (Otherwise just leave the column empty)
      vii. | Transfer the results for the salts into the original table:
         
           |   Add the BiGG ID of each ion to the component containing the ion such that in each row only one ion is represented.
           |   If duplicates occur only add the according ion to the first component containing it.

3. Check if the medium definition misses relevant components like water, iron, oxygen or carbondioxide.
4. If relevant components are missing, search for reasonable explanations to add the according components or justify why 
   these are missing.
5. The medium is now ready to be used for growth simulation! :) 
