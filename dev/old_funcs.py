# @TODO check for new-ish functionalities / merge with create_gp and delete
# @DEPRECATE
def create_gpr_from_locus_tag(model: libModel, locus_tag: str, email: str) -> tuple[GeneProduct, libModel]:
    """Creates GeneProduct in the given model

    **Deprecation warning**: will be deprecated in a future update.

    Args:
        - model (libModel): 
            Model loaded with libSBML
        - locus_tag (str): 
            NCBI compatible locus_tag
        - email (str): 
            User Email to access the NCBI Entrez database

    Returns:
        tuple: 
            libSBML GeneProduct (1) & libSBML model (2)

            (1) GeneProduct: Created gene product
            (2) libModel: Model containing the created gene product
    """
    mes = f'create_gpr_from_locus_tag will be deprecated in a future update.'
    warnings.warn(mes,type=FutureWarning)

    Entrez.email = email
    name, locus = search_ncbi_for_gpr(locus_tag)
    gpr = model.getPlugin(0).createGeneProduct()
    gpr.setName(name)
    gpr.setId(locus_tag)
    gpr.setMetaId('meta_' + locus_tag)
    gpr.setLabel(locus_tag)
    gpr.setSBOTerm("SBO:0000243")
    add_cv_term_genes(locus_tag, 'NCBI', gpr)
    return gpr, model

# @TODO : check, if the function (in ths way), is still used anywhere and adjust to the 
# new one
# @DEPRECATE
def old_create_gp(model: libModel, model_id: str, name: str, locus_tag: str, protein_id: str) -> tuple[GeneProduct, libModel]:
    """Creates GeneProduct in the given model

    Args:
        - model (libModel): 
            Model loaded with libSBML
        - model_id (str): 
            ID identical to ID that CarveMe adds from the NCBI FASTA input file
        - name (str): 
            Name of the GeneProduct
        - locus_tag (str): 
            Genome-specific locus tag used as label in the model
        - protein_id (str): 
            NCBI Protein/RefSeq ID

    Returns:
        tuple: 
            libSBML GeneProduct (1) & libSBML model (2)

            (1) GeneProduct: Created gene product
            (2) libModel: Model containing the created gene product
    """
    id_db = None
    gp = model.getPlugin(0).createGeneProduct()
    gp.setId(model_id) # libsbml advised to use set/getIdAttribute
    gp.setName(name)
    gp.setLabel(locus_tag)
    gp.setSBOTerm('SBO:0000243')
    gp.setMetaId(f'meta_{model_id}')
    if re.fullmatch(r'^(((AC|AP|NC|NG|NM|NP|NR|NT|NW|WP|XM|XP|XR|YP|ZP)_\d+)|(NZ_[A-Z]{2,4}\d+))(\.\d+)?$', protein_id, re.IGNORECASE):
        id_db = 'REFSEQ'
    elif re.fullmatch(r'^(\w+\d+(\.\d+)?)|(NP_\d+)$', protein_id, re.IGNORECASE): id_db = 'NCBI'
    if id_db: add_cv_term_genes(protein_id, id_db, gp)
    return gp, model


# @DISCUSSION: Add to metabolite handling for `.entities.build_metabolite_bigg`?
def add_charges_chemical_formulae_to_metabs(missing_metabs: pd.DataFrame) -> pd.DataFrame:
   """Adds charges & chemical formulae from CHEBI/BiGG to the provided dataframe

   Args:
      - missing_metabs (pd.DataFrame): 
         Table containing metabolites & the respective CHEBI & BiGG IDs
         
   Returns:
      pd.DataFrame: 
         Input table extended with the charges & chemical formulas obtained from CHEBI/BiGG
   """
   
   # Finds the charges through the ChEBI/BiGG API, defaults to: 0
   def find_charge(row: pd.Series) -> int:
      chebi_id = str(int(row.get('ChEBI'))) if not math.isnan(float(row.get('ChEBI'))) else None
      bigg_id = str(row.get('bigg_id'))
      charge = None
      if chebi_id:  # Get charge from ChEBI (Returns always a charge)
         chebi_entity = libchebipy.ChebiEntity('CHEBI:' + chebi_id)
         return chebi_entity.get_charge()
      elif bigg_id != 'nan':  # Get charge from BiGG if no ChEBI ID available
         try:
            charge = requests.get(BIGG_METABOLITES_URL + bigg_id[:-2]).json()['charges'][0]  # Take first charge
         except ValueError:
            pass   
         # If no charge was found, charge=0
         return charge if charge else 0
   
   # Finds the chemical formula through the ChEBI/BiGG API, defaults to: 'No formula'
   def find_formula(row: pd.Series) -> str:
      chebi_id = str(int(row.get('ChEBI'))) if not math.isnan(float(row.get('ChEBI'))) else None
      bigg_id, chem_form = str(row.get('bigg_id')), str(row.get('Chemical Formula'))
      chem_formula = None
      if chebi_id: # Get formula from ChEBI
         chebi_entity = libchebipy.ChebiEntity('CHEBI:' + chebi_id)
         chem_formula = chebi_entity.get_formula()
      if not chem_formula:  # If no formula was found with ChEBI/No ChEBI ID available
         if bigg_id != 'nan': # Get formula from BiGG
            try:
               chem_formula = requests.get(BIGG_METABOLITES_URL + bigg_id[:-2]).json()['formulae'][0]  # Take first formula
            except (ValueError, IndexError) as e:
               pass
         if not chem_formula: # If no formula was found with BiGG ID
            # Get formula already existing in dataframe or set to 'No formula'
            chem_formula = chem_form if chem_form != 'nan' else '*'
      return chem_formula
   
   missing_metabs['charge'] = missing_metabs.apply(find_charge, axis=1)
   missing_metabs['New Chemical Formula'] = missing_metabs.apply(find_formula, axis=1)
   missing_metabs['Chemical Formula'] = missing_metabs['New Chemical Formula']
   missing_metabs.drop('New Chemical Formula', axis=1, inplace=True)
   
   return missing_metabs


# @DISCUSSION
# @ASK Still required? Useful for something?
# -----------
def compare_bigg_ids(id1: str, id2: str) -> bool:
    """Compares two BiGG strings/IDs & Returns True if one BiGG ID matches most of the other

    Args:
        - id1 (str): 
            ID 1
        - id2 (str): 
            ID 2

    Returns:
        bool: 
            Indicates if most of one string contained in the other
    """
    id1_split, id2_split, id1_single_comp, id2_single_comp, id1_comp, id2_comp = None, None, None, None, None, None
    
    if '_' in id1: id1_split = re.split(r'_([a-zA-Z]|[0-9])$', id1)[0]
    if '_' in id2: id2_split = re.split(r'_([a-zA-Z]|[0-9])$', id2)[0]
    if id1.endswith(ALL_BIGG_COMPARTMENTS_ONE_LETTER): id1_single_comp = id1[:-1]
    if id2.endswith(ALL_BIGG_COMPARTMENTS_ONE_LETTER): id2_single_comp = id2[:-1]
    if id1.endswith(ALL_BIGG_COMPARTMENTS_TWO_LETTER): id1_comp = id1[:-2]
    if id2.endswith(ALL_BIGG_COMPARTMENTS_TWO_LETTER): id2_comp = id2[:-2]
    
    similar_ids = False
    if id1 == id2: similar_ids = True  # Both IDs are same
    
    elif id1_split and id2_split and (id1_split == id2_split): similar_ids = True # Both IDs are same but from different compartments
    elif id2_split and (id1 == id2_split): similar_ids = True # - "" -
    elif id1_split and (id1_split == id2): similar_ids = True # - "" -
    
    elif id1_single_comp and id2_single_comp and (id1_single_comp == id2_single_comp): similar_ids = True
    elif id2_single_comp and (id1 == id2_single_comp): similar_ids = True
    elif id1_single_comp and (id1_single_comp == id2_single_comp): similar_ids = True 
    
    elif id1_comp and id2_comp and (id1_comp == id2_comp): similar_ids = True
    elif id2_comp and (id1 == id2_comp): similar_ids = True
    elif id1_comp and (id1_comp == id2): similar_ids = True
    
    elif id1_split and id2_single_comp and (id1_split == id2_single_comp): similar_ids = True
    elif id2_split and id1_single_comp and (id1_single_comp == id2_split): similar_ids = True
    
    elif id1_split and id2_comp and (id1_split == id2_comp): similar_ids = True
    elif id2_split and id1_comp and (id1_comp == id2_split): similar_ids = True
    
    elif id1_comp and id2_single_comp and (id1_comp == id2_single_comp): similar_ids = True
    elif id2_comp and id1_single_comp and (id1_single_comp == id2_comp): similar_ids = True
    
    else: similar_ids = False

    return similar_ids

# @TEST
# @NOTE : A lot of warnings
# @DEPRECATE: Nice function but not used anymore & will have issues running
def keep_only_bigg_reactions_in_certain_compartments(complete_df: pd.DataFrame) -> pd.DataFrame:
    """Extracts all possible BiGG ID variations from database for a BiGG reaction ID, gets the metabolite compartments 
    & returns table containing only reactions which happen in one of the provided compartments
        
    Args:
        - complete_df (pd.DataFrame): 
            Table containing at least the column 'bigg_id'.
        
    Returns:
        pd.DataFrame: 
            Table containing reactions & their compartments
    """
    tqdm.pandas()
    
    # (1) Find all occurrencs of a BiGG reaction ID in bigg_reactions table in database
    def get_all_similar_bigg_ids(bigg_id_in: str) -> list[str]:
        
        if '_' in bigg_id_in: bigg_id = re.split(r'_([a-zA-Z]|[0-9])$', bigg_id_in)[0]
        elif bigg_id_in.endswith(ALL_BIGG_COMPARTMENTS_ONE_LETTER): bigg_id = bigg_id_in[:-1]
        elif bigg_id_in.endswith(ALL_BIGG_COMPARTMENTS_TWO_LETTER): bigg_id = bigg_id_in[:-2]
        else: bigg_id = bigg_id_in
        
        query = f"SELECT id, INSTR(id, '{bigg_id}') bi FROM bigg_reactions WHERE bi > 0"
        result = con.execute(query).fetchall()
        result = [result_tuple[0] for result_tuple in result] if result else [bigg_id_in]
        result = [res for res in result if compare_bigg_ids(bigg_id, res)]
        return result
    
    # (2) Use list of all BiGG IDs obtained from database table bigg_reactions to get 'metabolites'
    # get_react_compartment moved to outer scope due to multiprocessing Pool (see line 95)
    def multi_get_reaction_compartment(complete_df: pd.DataFrame) -> list:
        """Takes a dataframe and runs get_reaction_compartment() in multiple
        processes on the 'bigg_id' column

        Args:
            - complete_df (pd.DataFrame): Table containing at least the columns 'bigg_id' & 'KEGG'/'BioCyc'

        Returns:
            list: List of compartments
        """
        with Pool() as pool:
            results = []
            for out in tqdm(pool.imap(get_reaction_compartment_from_bigg, complete_df.loc[:, "bigg_id"], chunksize=20), total=len(complete_df)):
                results.append(out)

        return results

    print(complete_df.columns)
    
    # Connect to database & get similar IDs (1)
    print('Getting all similar IDs...')
    con = sqlite3.connect(PATH_TO_DB)  # Open connection to database
    complete_df.loc[:, 'bigg_id_list'] = complete_df.loc[:, 'bigg_id'].progress_map(get_all_similar_bigg_ids)
    # complete_df.progress_apply(get_all_similar_bigg_ids, axis=1)
    con.close()  # Close connection to database

    # Adjust table to contain one BiGG ID per row from bigg_id_list (1)
    complete_df.loc[:, 'id_group'] = complete_df['bigg_id'].ne(complete_df['bigg_id'].shift()).cumsum()  # Group similar IDs
    complete_df.drop(labels='bigg_id', axis=1, inplace=True)  # Drop 'bigg_id' as no longer required
    complete_df = complete_df.explode('bigg_id_list', ignore_index=True)  # Expand 'bigg_id_list' column
    complete_df.rename(columns={'bigg_id_list': 'bigg_id'}, inplace=True)  # Rename 'bigg_id_list' to 'bigg_id'

    # (2) Get all compartments for each reaction from BiGG database API
    print(f'Getting all IDs with correct compartment {VALID_COMPARTMENTS.keys()}...')
    results = multi_get_reaction_compartment(complete_df)
    complete_df["compartment"] = results

    # complete_df.progress_apply(get_reaction_compartment, axis=1)  # (2)

    # (3) Remove reactions with compartment = NaN
    complete_df.dropna(subset=['compartment'], inplace=True)

    return complete_df

# @TEST
# @DEPRECATE: Nice function but not used anymore & might have issues running
def get_bigg_db_mapping(map_to:str='BioCyc', metabolites:bool=True) -> pd.DataFrame:
    """Download a mapping of BiGG IDs to a specified database.

    Args:
        map_to (str, optional): 
            Name of the database to map to. 
            Ideally a column of the table in the database, 
            but SEED, KEGG and BioCyc are valid as well. 
            Defaults to 'BioCyc'.
        metabolites (bool, optional): 
            Flag to map reaction (False) or metabolite (True) IDs. 
            Defaults to True.

    Raises:
        - KeyError: Given database name not found in database. Cannot perform mapping.

    Returns:
        pd.DataFrame: 
            The mapping as a table.
    """

    # adjust name to map to if necessary
    reac_or_comp = 'Compound' if metabolites else 'Reaction'
    table_name = 'bigg_metabolites' if metabolites else 'bigg_reactions'
    if map_to in ['SEED','KEGG','Reactome']:
        map_to =  ' '.join([map_to, reac_or_comp])
    
    # download BiGG tables from database
    # ----------------------------------
    # build connection to DB
    connection = sqlite3.connect(PATH_TO_DB)
    cursor = connection.cursor()

    # retrieve only mappings to a specific database
    result = cursor.execute('SELECT 1 FROM PRAGMA_TABLE_INFO(?) WHERE name = ?',(table_name,map_to))
    possible_db = result.fetchone()
    if possible_db:
        query = f'SELECT * FROM {table_name} WHERE {map_to} IS NOT NULL'
    else:
        raise KeyError('Given database name not found in database. Cannot perform mapping.')

    # actually load data
    data = load_a_table_from_database(query)
    data = data.explode(map_to, ignore_index=True)

    # reduce columns to mapping only
    data = data[['id',map_to]]
    data.rename(columns={'id':'bigg_id'}, inplace=True)

    # filter for compartment in case of reactions
    if not metabolites:
        data = keep_only_bigg_reactions_in_certain_compartments(data)

    # close connection to database
    connection.close()
    
    return data
