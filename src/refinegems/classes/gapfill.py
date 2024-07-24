"""_summary_
"""

__author__ = "Famke Baeuerle, Gwendolyn O. Döbel, Carolin Brune and Dr. Reihaneh Mostolizadeh"

############################################################################
# requirements
############################################################################

from abc import ABC, abstractmethod

import pandas as pd
import numpy as np

from libsbml import Model as libModel
from bioservices.kegg import KEGG
import io
import re
from tqdm import tqdm
tqdm.pandas()

from ..curation.db_access.db import get_bigg_db_mapping, compare_bigg_model, add_stoichiometric_values_to_reacs
from ..utility.io import load_a_table_from_database
from ..utility.entities import get_model_reacs_or_metabs
from ..curation.db_access.kegg import parse_KEGG_gene, parse_KEGG_ec

############################################################################
# variables
############################################################################

############################################################################
# classes
############################################################################

# -------------
# abtract class
# -------------

class GapFiller(ABC):

    def __init__(self) -> None:
        self.full_gene_list = None
        self.geneid_type = 'ncbi' # @TODO
        self._statistics = dict()
        
    #@abstractmethod
    #def load_gene_list(self):
    #    pass

    @abstractmethod
    def get_missing_genes(self, model):
        pass
    
    @abstractmethod
    def get_missing_reacs(self,model):
        pass
    
    # @abstractmethod
    # def run(self,model):
    #     pass

    def fill_model(self, model):
        pass

    def calculate_stats(self):
        pass

    def report(self):
        pass

# --------------------
# Gapfilling with KEGG
# --------------------

class KEGGapFiller(GapFiller):

    def __init__(self, organismid) -> None:
        super().__init__()
        self.organismid = organismid
        # self.report = dict()
        self.manual_curation = dict()
        
    # @TODO: progress bar and parallelising
    # @TODO: logging
    def get_missing_genes(self, model:libModel):
    
        # Function originally from refineGEMs.genecomp/refineGEMs.KEGG_analysis/entities --- Modified
        def get_model_genes(model: libModel) -> pd.DataFrame:
            """Extracts KEGG Genes from given model

            Args:
                - model (model-libsbml): 
                    Model loaded with libSBML

            Returns:
                pd.DataFrame: 
                    Table with all KEGG Genes in the model
            """
            genes_in_model = []
            for gene in model.getPlugin(0).getListOfGeneProducts():
                cv_terms = gene.getCVTerms()
                if cv_terms:
                    for cv_term in cv_terms:
                        for idx in range(cv_term.getNumResources()):
                            uri = cv_term.getResourceURI(idx)
                            if 'kegg.genes' in uri: 
                                genes_in_model.append(re.split('kegg.genes:|kegg.genes/',uri)[1]) # work with olf/new pattern

            return pd.DataFrame(genes_in_model, columns=['orgid:locus'])
        
        
        # Step 1: get genes from model
        # ----------------------------
        genes_in_model = get_model_genes(model)
        
        # Step 2: get genes of organism from KEGG
        # ---------------------------------------
        gene_KEGG_list = KEGG().list(self.organismid)
        gene_KEGG_table = pd.read_table(io.StringIO(gene_KEGG_list), header=None)
        gene_KEGG_table.columns = ['orgid:locus','CDS','position','protein']
        gene_KEGG_table = gene_KEGG_table[['orgid:locus']]
        
        # Step 3: KEGG vs. model genes -> get missing genes for model
        # ----------------------------
        genes_not_in_model = gene_KEGG_table[~gene_KEGG_table['orgid:locus'].isin(genes_in_model['orgid:locus'])]
        
        # Step 4: extract locus tag
        # -------------------------
        genes_not_in_model['locus_tag'] = genes_not_in_model['orgid:locus'].str.split(':').str[1]
        
        # Step 5: map to EC via KEGG
        # --------------------------
        # @DEBUGGING ...................
        # genes_not_in_model = genes_not_in_model.iloc[330:350,:]
        # ..............................
        geneKEGG_mapping = pd.DataFrame.from_dict(list(genes_not_in_model['orgid:locus'].progress_apply(parse_KEGG_gene)))
        genes_not_in_model = genes_not_in_model.merge(geneKEGG_mapping, how='left', on='orgid:locus')
        
        # @TODO : What to report where and when
        # self.report['missing genes (total)'] = len(genes_not_in_model)
        
        return genes_not_in_model 
    
    # @TODO : logging
    # @TODO : paralellising possibilities?
    # @TODO : progress bar
    def get_missing_reacs(self, model, genes_not_in_model):
        
        # Step 1: get reactions from model 
        # --------------------------------
        reac_model_list = model.getListOfReactions()
        reac_model_list = [_.id[2:] for _ in reac_model_list] # crop 'R_' prefix

        # Step 2: filter missing gene list + extract ECs
        # ----------------------------------------------
        reac_options = genes_not_in_model[['ec-code','ncbiprotein']]        # get relevant infos for reacs
        missing_reacs = reac_options[['ec-code','ncbiprotein']].dropna()    # drop nas
        self.manual_curation['genes'] = reac_options.loc[~reac_options.index.isin(missing_reacs.index)]
        # check, if any automatic gapfilling is possible
        if len(missing_reacs) == 0:
            return None
        # transform table into EC-number vs. list of NCBI protein IDs
        eccode = missing_reacs['ec-code'].apply(pd.Series).reset_index().melt(id_vars='index').dropna()[['index', 'value']].set_index('index')
        ncbiprot = missing_reacs['ncbiprotein'].apply(pd.Series).reset_index().melt(id_vars='index').dropna()[['index', 'value']].set_index('index')
        missing_reacs = pd.merge(eccode,ncbiprot,left_index=True, right_index=True).rename(columns={'value_x':'ec-code','value_y':'ncbiprotein'})
        missing_reacs = missing_reacs.groupby(missing_reacs['ec-code']).aggregate({'ncbiprotein':'unique'}).reset_index()
        
        # Step 3: Map to MetaNetX
        # -----------------------
        # convert table into one EC-number a row
        mnx_reac_prop = load_a_table_from_database('mnx_reac_prop',False)
        mnx_reac_prop.drop('is_balanced', inplace=True, axis=1)
        mnx_reac_prop['ec-code'] = mnx_reac_prop['ec-code'].apply(lambda x: x.split(';') if isinstance(x,str) else None)
        mnx_reac_prop = mnx_reac_prop.explode('ec-code').dropna(subset='ec-code')
        # merge tables on EC-number
        reacs_mapped = missing_reacs.merge(mnx_reac_prop, on='ec-code', how='left')
        reacs_mapped.rename({'mnx_equation':'equation'}, inplace=True, axis=1)
        reacs_mapped['via'] = reacs_mapped['id'].apply(lambda x: 'MetaNetX' if x else None)
                
        # Step 4: map to BiGG
        # -------------------
        if reacs_mapped.id.isnull().values.any():
            # load BiGG reaction namespace
            bigg_reacs = load_a_table_from_database('bigg_reactions',False)
            bigg_reacs.dropna(subset='EC Number', inplace=True)
            bigg_reacs = bigg_reacs[['id','reaction_string','EC Number']].rename({'reaction_string':'equation','EC Number':'ec-code'}, inplace=False, axis=1)
            bigg_reacs['ec-code'] = bigg_reacs['ec-code'].apply(lambda x: x.split(',') if isinstance(x,str) else None)
            bigg_reacs = bigg_reacs.explode('ec-code')
            # merge unmapped entries with BiGG on EC-number 
            bigg_mapping = reacs_mapped[reacs_mapped['id'].isnull()][['ec-code','ncbiprotein']].merge(bigg_reacs, on=['ec-code'], how='left')
            bigg_mapping['via'] = bigg_mapping['id'].apply(lambda x: 'BiGG' if x else None)
            # combine with MNX results
            reacs_mapped = pd.concat([reacs_mapped[~reacs_mapped['id'].isnull()],bigg_mapping], axis=0, ignore_index=True)

        # Step 5: Map to KEGG
        # ------------------- 
        if reacs_mapped.id.isnull().values.any():
            # get KEGG EC number information
            kegg_mapped = pd.DataFrame.from_dict(list(reacs_mapped[reacs_mapped['id'].isnull()]['ec-code'].progress_apply(parse_KEGG_ec)))
            # kegg_mapped = kegg_mapped.explode('id')      
            kegg_mapped['is_transport'] = None
            kegg_mapped['via'] = kegg_mapped['id'].apply(lambda x: 'KEGG' if x else None)
            # join with NCBI protein ID
            kegg_mapped = reacs_mapped[reacs_mapped['id'].isnull()][['ec-code','ncbiprotein']].merge(kegg_mapped, on=['ec-code'], how='left')
            # merge with input
            reacs_mapped = pd.concat([reacs_mapped[~reacs_mapped['id'].isnull()],kegg_mapped], axis=0, ignore_index=True)
            reacs_mapped.mask(reacs_mapped.isna(), other=None, inplace=True)
                
        # Results
        # -------
        reacs_mapped = reacs_mapped.explode('id', ignore_index=True).explode('id', ignore_index=True)
        # need manual curation
        self.manual_curation['reacs'] = reacs_mapped[reacs_mapped['id'].isnull()]
                
        return reacs_mapped[~reacs_mapped['id'].isnull()] 
    
    
# ----------------------
# Gapfilling with BioCyc
# ----------------------

class BioCycGapFiller(GapFiller):
    """
    | @TODO: Write correct doc string
    |
    |
    |
    |These three TXT files can be obtained through creating SmartTables in BioCyc and exporting these as 
    |Spreadsheets with the parameter FrameID (one & two) or Common name (three). 
    |SmartTable one: Should contain 'Accession-2' and 'Reactions of gene', 
    |SmartTable two: Should contain 'Reaction', 'Reactants of reaction', 'Products of reaction', 'EC-Number', 'Reaction-Direction' and 'Spontaneous?'
    |SmartTable three: Should contain 'Compound', 'Chemical Formula' and 'InChI-Key'.

    Args:
        - GapFiller (_type_): 
            _description_
    """

    def __init__(self, model: libModel, biocyc_gene_tbl_path: str, biocyc_reacs_tbl_path: str) -> None:
        super().__init__()
        self._biocyc_gene_tbl = self.biocyc_gene_tbl(biocyc_gene_tbl_path) # @TODO: Das ist eigtl self.full_gene_list von GapFiller!
        self._biocyc_rxn_tbl = self.biocyc_rxn_tbl(biocyc_reacs_tbl_path)
        self._model = model
        self.missing_biocyc_genes = None
        self.missing_biocyc_reacs = None
        self.missing_biocyc_metabs = None

    @property
    def biocyc_gene_tbl(self):
        return self._BioCyc_gene_tbl
    
    @biocyc_gene_tbl.setter
    # @TODO: Hier sollten wir noch diskutieren, ob Accession-2 oder Accession-1 hard coded sein sollten oder nicht
    #        Dafür müssten wir nochmal ein paar entries in BioCyc dazu anschauen.
    # Locus tags in GenBank GFF file == BioCyc Accession-2 == Old locus tags in RefSeq GFF file
    # Locus tags in RefSeq GFF file == BioCyc Accession-1
    # Label in model == Locus tag from GenBank GFF file == BioCyc Accession-2
    def biocyc_gene_tbl(self, biocyc_gene_tbl_path: str) -> pd.DataFrame:
        """Parses TSV file from BioCyc to retrieve 'Accession-2' & the corresponding 'Reactions of gene'
    
        Args:
           - inpath (str): 
              Path to file from BioCyc containing the 'Accession-2' to 'Reactions of gene' mapping
           
        Returns:
           pd.DataFrame: 
              Table containing only rows where a 'Reaction of gene' exists
        """
        biocyc_genes = pd.read_table(biocyc_gene_tbl_path, usecols=['Accession-2', 'Reactions of gene'], dtype=str)
        biocyc_genes.rename(columns={'Accession-2': 'locus_tag', 'Reactions of gene': 'Reaction'}, inplace=True)
        biocyc_genes.replace('', np.nan, inplace=True)
        biocyc_genes.dropna(inplace=True)
        self._BioCyc_gene_tbl = biocyc_genes

    @property
    def biocyc_rxn_tbl(self):
        return self._biocyc_rxn_tbl

    @biocyc_rxn_tbl.setter
    def biocyc_rxn_tbl(biocyc_reacs_tbl_path: str) -> pd.DataFrame:
        """Parses TSV file from BioCyc to retrieve 'Reaction', 'Reactants of reaction', 'Products of reaction', 'EC-Number',
           'Reaction-Direction' & 'Spontaneous?'

        Args:
           - biocyc_reacs_tbl_path (str):   
              Path to file from BioCyc containing the following columns:
              'Reaction' 'Reactants of reaction' 'Products of reaction' 'EC-Number' 'Reaction-Direction' 
              'Spontaneous?'

        Returns:
           pd.DataFrame: 
              Table containing all BioCyc reactions from provided file
        """
   
        biocyc_reacs = pd.read_table(biocyc_reacs_tbl_path, usecols=
                                     ['Reaction', 'Reactants of reaction', 'Products of reaction', 'EC-Number', 
                                      'Reaction-Direction', 'Spontaneous?'],
                                     dtype=str
                                     )
        biocyc_reacs.rename(columns=
                            {'Reactants of reaction': 'Reactants', 'Products of reaction': 'Products', 'EC-Number': 'ec-code'},
                            inplace=True
                            )
        biocyc_reacs.replace('', np.nan, inplace=True)
        biocyc_reacs['Spontaneous?'] = biocyc_reacs['Spontaneous?'].fillna('F')
        biocyc_reacs.dropna(subset=['Reaction', 'Reactants', 'Products', 'EC'], inplace=True)
        return biocyc_reacs
    
    def get_missing_genes(self):
        """Retrieves the missing genes and reactions from the BioCyc table according to the 'Accession-2' identifiers
        """
        
        # Step 1: get genes from model
        # ----------------------------
        geneps_in_model = [_.getLabel() for _ in self._model.getPlugin(0).getListOfGeneProducts()]

        # Step 2: Get genes of organism from BioCyc
        # -----------------------------------------
        # For now see setter: BioCyc_gene_tbl
        
        # Step 3: BioCyc vs. model genes -> get missing genes for model
        # -------------------------------------------------------------
        self.missing_biocyc_genes = self.biocyc_gene_tbl[~self.biocyc_gene_tbl['locus_tag'].isin(geneps_in_model['locus_tag'])]

        # Step 4: Get amount of missing genes from BioCyc for statistics
        # --------------------------------------------------------------
        self._statistics['Protein']['Total'] = len(self.missing_biocyc_genes['locus_tag'].unique().tolist())

    # @TODO Result with columns: ec-code | ncbiprotein | id | equation (Gibt es das in BioCyc?) (| reference )| is_transport (Gibt es das in BioCyc?) | via
    def get_missing_reacs(self) -> tuple[tuple[pd.DataFrame, pd.DataFrame], pd.DataFrame]:
        """Subsets the BioCyc table with the following columns:
         'Reaction' 'Reactants of reaction' 'Products of reaction' 'EC-Number' 'Reaction-Direction' 'Spontaneous?'
         to obtain the missing reactions with all the corresponding data
        """
        
        # Step 1: get reactions from model
        # --------------------------------
        reac_model_list = get_model_reacs_or_metabs(self._model)

        # Step 2: filter missing gene list + extract ECs
        # ----------------------------------------------
        # Expand missing genes result table to merge with Biocyc reactions table
        missing_biocyc_genes = pd.DataFrame(
           self.missing_biocyc_genes['Reaction'].str.split('//').tolist(), index=self.missing_biocyc_genes['locus_tag']
           ).stack()
        missing_biocyc_genes = missing_biocyc_genes.reset_index([0, 'locus_tag'])
        missing_biocyc_genes.columns = ['locus_tag', 'Reaction']
        missing_biocyc_genes['Reaction'] = missing_biocyc_genes['Reaction'].str.strip()

        # Get missing reactions from missing genes
        missing_reacs = self.missing_biocyc_genes.merge(self.biocyc_rxn_tbl, on='Reaction')

        # Turn entries with '//' into lists
        missing_reacs['Reactants'] = missing_reacs['Reactants'].str.split('\s*//\s*')
        missing_reacs['Products'] = missing_reacs['Products'].str.split('\s*//\s*')
        missing_reacs['ec-code'] = missing_reacs['EC'].str.split('\s*//\s*')

        # Step 3: Get content for column ncbiprotein
        # ------------------------------------------
        # Turn locus_tag column into lists of locus tags per reaction
        locus_tags_as_list = missing_reacs.groupby('Reaction')['locus_tag'].apply(list).reset_index(name='locus_tag')
        missing_reacs.drop('locus_tag', axis=1, inplace=True)
        missing_reacs = locus_tags_as_list.merge(missing_reacs, on='Reaction')

        # @TODO: Map locus tags to corresponding protein IDs

        # Step 4: Get amount of missing reactions from BioCyc for statistics
        # ------------------------------------------------------------------
        self._statistics['Reaction']['Total'] = len(missing_reacs['Reaction'].unique().tolist())

        
        return missing_reacs

    
# ----------------
# Gapfilling no DB
# ----------------
# @NOTE: Ideas from Gwendolyn O. Döbel for lab strains
# Get all  possible genes by filtering .gff according to 'bio_type=protein_coding' & 'product=hypothetical protein'
# Compare the list of genes with the ones already in the model & add all missing genes
# Before adding to model check if for all genes that are missing for IMITSC147 identifiers exist
# -> Create tables mapping locus tag to old ID, locus tag to new ID & merge 
# -> Specify user input locus_tag start from NCBI PGAP

class GeneGapFiller(GapFiller):
    pass

############################################################################
# functions
############################################################################


############################################################################
# For filtering
############################################################################

# Get BiGG BioCyc
bigg2biocyc_reacs = get_bigg_db_mapping('BioCyc',False)

# Subset missing_reactions with BiGG BioCyc
missing_reactions.rename(columns={'Reaction': 'BioCyc'}, inplace=True)
# TODO: collect non-matched entries // return all // namespace independance ????
missing_reactions = bigg2biocyc_reacs.merge(missing_reactions, on='BioCyc')

# Get amount of missing reactions that have a BiGG ID
self._statistics['Reaction']['Have BiGG ID'] = len(missing_reactions['BioCyc'].unique().tolist())

# Subset missing_reactions with model_reacs
missing_reactions = compare_bigg_model(missing_reactions, model_reacs)

# Get amount of missing reactions that are not in the model
self._statistics['Reaction']['Can be added'] = len(missing_reactions['bigg_id'].unique().tolist())

# Add reactants & products dictionary with stoichiometric values to the reactions table
missing_reactions = add_stoichiometric_values_to_reacs(missing_reactions)

# Get all metabolites for the missing reactions
biocyc_metabs_from_reacs, bigg_metabs_from_reacs = extract_metabolites_from_reactions(missing_reactions)
