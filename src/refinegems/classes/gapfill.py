"""_summary_
"""

__author__ = "Famke Baeuerle, Gwendolyn O. Döbel, Carolin Brune and Dr. Reihaneh Mostolizadeh"

############################################################################
# requirements
############################################################################

from abc import ABC, abstractmethod

import cobra
import io
import numpy as np
import pandas as pd
import re
import warnings

from bioservices.kegg import KEGG
from libsbml import Model as libModel
from typing import Literal, Union

from tqdm import tqdm
tqdm.pandas()

from ..curation.db_access.db import get_ec_from_ncbi, get_ec_via_swissprot
from ..utility.io import load_a_table_from_database, parse_fasta_headers, parse_gff_for_cds
from ..utility.entities import get_model_reacs_or_metabs, create_gp, create_gpr
from ..curation.db_access.kegg import parse_KEGG_gene, parse_KEGG_ec

############################################################################
# variables
############################################################################


############################################################################
# where to put 
############################################################################

# Mapping of EC numbers
# ---------------------

def map_ec_to_reac(table, use_MNX=True, use_BiGG=True, use_KEGG=True):
    
    def _map_ec_to_reac_mnx(unmapped_reacs):
    
        # input: pd.DataFrame with a least the ec-code column
        # load MNX reac prop table
        mnx_reac_prop = load_a_table_from_database('mnx_reac_prop',False)
        # convert table into one EC-number a row
        mnx_reac_prop.drop('is_balanced', inplace=True, axis=1)
        mnx_reac_prop['ec-code'] = mnx_reac_prop['ec-code'].apply(lambda x: x.split(';') if isinstance(x,str) else None)
        # exclude entries without EC-number
        mnx_reac_prop = mnx_reac_prop.explode('ec-code').dropna(subset='ec-code')
        # merge with unmapped reactions
        reacs_mapped = unmapped_reacs.merge(mnx_reac_prop, on='ec-code', how='left')
        reacs_mapped.rename({'mnx_equation':'equation'}, inplace=True, axis=1)
        reacs_mapped['via'] = reacs_mapped['id'].apply(lambda x: 'MetaNetX' if x else None)
        
        return reacs_mapped      
    

    def _map_ec_to_reac_bigg(unmapped_reacs):
        
        # load BiGG reaction namespace
        bigg_reacs = load_a_table_from_database('bigg_reactions',False)
        bigg_reacs.dropna(subset='EC Number', inplace=True)
        bigg_reacs = bigg_reacs[['id','reaction_string','EC Number']].rename({'reaction_string':'equation','EC Number':'ec-code'}, inplace=False, axis=1)
        bigg_reacs['ec-code'] = bigg_reacs['ec-code'].apply(lambda x: x.split(',') if isinstance(x,str) else None)
        bigg_reacs = bigg_reacs.explode('ec-code')
        # merge with unmapped reactions
        bigg_mapping = unmapped_reacs.merge(bigg_reacs, on=['ec-code'], how='left')
        bigg_mapping.mask(bigg_mapping.isna(), other=None, inplace=True)
        # make conform to format
        bigg_mapping['reference'] = None
        bigg_mapping['is_transport'] = None
        bigg_mapping['via'] = bigg_mapping['id'].apply(lambda x: 'BiGG' if x else None)
        
        return bigg_mapping


    def _map_ec_to_reac_kegg(unmapped_reacs):
        
        # get KEGG EC number information
        kegg_mapped = pd.DataFrame.from_dict(list(unmapped_reacs.progress_apply(parse_KEGG_ec)))   
        kegg_mapped['is_transport'] = None
        kegg_mapped['via'] = kegg_mapped['id'].apply(lambda x: 'KEGG' if x else None)
        kegg_mapped.explode(column='id')
        
        return kegg_mapped


    # input table should have format: 
    #   ec-code | ncbiprotein
    #   one EC number per row, list of ncbiprotein per row allowed
    if len(table.columns) != 2 or 'ec-code' not in table.columns or 'ncbiprotein' not in table.columns:
        raise ValueError('Wrong table format. Cannot map EC to reaction.')
    
    # map to MetaNetX
    if use_MNX:
        table = _map_ec_to_reac_mnx(table)
    
    # map to BiGG
    if use_BiGG:
        if 'id' in table.columns:
            to_map = table[table['id'].isna()][['ec-code','ncbiprotein']]    
            if len(to_map) > 0:        
                to_map = _map_ec_to_reac_bigg(to_map) 
                table = pd.concat([to_map,table[~table['id'].isna()]])
        else:
            table = _map_ec_to_reac_bigg(table)
    
    # map to KEGG
    if use_KEGG: 
        if 'id' in table.columns:
            to_map = table[table['id'].isna()][['ec-code','ncbiprotein']] 
            if len(to_map) > 0:
                to_map = _map_ec_to_reac_kegg(to_map) 
                table = pd.concat([to_map,table[~table['id'].isna()]])
        else:
            table = _map_ec_to_reac_kegg(table)
    
    # explode
    table = table.explode('id', ignore_index=True).explode('id', ignore_index=True)
        
    # output:   ec-code ncbiprotein	id	equation	reference	is_transport	via
    return table 


############################################################################
# classes
############################################################################

# -------------
# abtract class
# -------------

class GapFiller(ABC):

    def __init__(self) -> None:
        self.full_gene_list = None  # really a good idea? GeneGapFiller does not need it at all
        self.geneid_type = 'ncbi' # @TODO
        self._statistics = dict()
        self.manual_curation = dict()
        
    #@abstractmethod
    #def load_gene_list(self):
    #    pass
    
    # abstract methods
    # ----------------

    @abstractmethod
    def get_missing_genes(self, model):
        pass
    
    @abstractmethod
    def get_missing_reacs(self,model):
        pass
    
    # finding the gaps
    # ----------------

    def _find_reac_in_model(self, model: cobra.Model, eccode:str, id:str, 
               idtype:Literal['MetaNetX','KEGG','BiGG', 'BioCyc'], 
               include_ec_match:bool=False) -> Union[None, list]:
        # @TODO Ensure that user has requested BioCyc identifiers.org version? 
        # -> Could be done with polish_annotations
        MAPPING = {
            'MetaNetX':'metanetx.reaction', 
            'KEGG':'kegg.reaction',
            'BiGG':'bigg.reaction',
            'BioCyc': 'metacyc.reaction'
            }

        found = []
        for r in model.reactions:
            if MAPPING[idtype] in r.annotation.keys():
                if (isinstance(r.annotation[MAPPING[idtype]],list) and 
                    id in r.annotation[MAPPING[idtype]]):
                    found.append(r.id)
                elif (isinstance(r.annotation[MAPPING[idtype]],str) and 
                      id == r.annotation[MAPPING[idtype]]):
                    found.append(r.id)
            if include_ec_match and eccode and 'ec-code' in r.annotation.keys():
                if (isinstance(r.annotation['ec-code'],list) and 
                    eccode in r.annotation['ec-code']):
                    found.append(r.id)
                elif (isinstance(r.annotation['ec-code'],str) and 
                      eccode == r.annotation['ec-code']):
                    found.append(r.id)
            
        found = list(set(found))
        
        if len(found) > 0:
            return found
        
        return None
    
    
    # actual "Filling" part
    # ---------------------
    
    # @TODO logging? or remove self?
    def add_genes_from_table(self,model:libModel, gene_table:pd.DataFrame) -> None:
        """Create new GeneProduct for a table of genes in the format:
        
        | ncbiprotein | locus_tag | UniProt | ... |
        
        The dots symbolise additional columns, that can be passed to the function, 
        but will not be used by it. The other columns, except UniProt, are required.

        Args:
            - model (libModel): 
                _description_
            - gene_table (pd.DataFrame): 
                The table with the genes to add. At least needs the columns
                '','' and 'ec-code'.
        """
    
        # ncbiprotein | locus_tag | ...
        # work on a copy to ensure input stays the same
        gene_table = gene_table.copy()
        # gene_table.drop(columns=['ec-code'],inplace=True)
        
        # create gps from the table and add them to the model
        if 'UniProt' in gene_table.columns:
            for idx,x in gene_table.iterrows():
                create_gp(model, x['ncbiprotein'], 
                        locus_tag=x['locus_tag'],
                        uniprot=(x['UniProt'],True))
        else:
            for idx,x in gene_table.iterrows():
                create_gp(model, x['ncbiprotein'], 
                        locus_tag=x['locus_tag'])
                

    # @TODO seems very ridgid, better ways to find the ids?
    def add_gene_reac_associations_from_table(model:libModel,
                                            reac_table:pd.DataFrame) -> None:
        """Using a table with at least the columns 'ncbiprotein' 
        (containing e.g. NCBI protein identifier (lists), should be gene IDs in the model)
        and 'add_to_GPR' (containing reactions identifier (lists)), add the gene IDs to the 
        GPRs of the corresponding reactions.

        Args:
            - model (libModel): 
                The model loaded with libSBML.
            - reac_table (pd.DataFrame): 
                The table containing at least the columns 'ncbiprotein' (gene IDs) and
                'add_to_GPR' (reaction IDs)
        """
        
        model_gene_ids = [_.getId() for _ in model.getPlugin(0).getListOfGeneProducts()]
        
        # get each unique ncbiprotein vs reaction mapping
        reac_table = reac_table[['ncbiprotein','add_to_GPR']]
        reac_table = reac_table.explode('ncbiprotein').explode('add_to_GPR')
        reac_table.drop_duplicates(inplace=True)
        
        # add the genes to the corresponding GPRs
        for idx,row in reac_table.iterrows():
            # check, if G_+ncbiprotein in model
            # if yes, add gpr
            geneid = 'G_'+row['ncbiprotein'].replace('.','_')
            reacid = 'R_'+row['add_to_GPR']
            if geneid in model_gene_ids:
                create_gpr(model.getReaction(reacid),geneid)
            # else, print warning
            else:
                mes = f'Cannot find {geneid} in model. Should be added to {reacid}'
                warnings.warn(mes,UserWarning)
    



    def fill_model(self, model):
        pass
    
    # reporting
    # ---------

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
        # print(UserWarning('Running in debugging mode.'))
        # ..............................
        geneKEGG_mapping = pd.DataFrame.from_dict(list(genes_not_in_model['orgid:locus'].progress_apply(parse_KEGG_gene)))
        genes_not_in_model = genes_not_in_model.merge(geneKEGG_mapping, how='left', on='orgid:locus')
        
        # @TODO : What to report where and when
        # self.report['missing genes (total)'] = len(genes_not_in_model)
        
        return genes_not_in_model 
    
    # @TODO : logging
    # @TODO : paralellising possibilities?
    # @TODO : progress bar
    def get_missing_reacs(self,model:cobra.Model,genes_not_in_model):
 
        # Step 1: filter missing gene list + extract ECs
        # ----------------------------------------------
        reac_options = genes_not_in_model[['ec-code','ncbiprotein']]        # get relevant infos for reacs
        missing_reacs = reac_options[['ec-code','ncbiprotein']].dropna()    # drop nas
        # self.manual_curation['genes'] = reac_options.loc[~reac_options.index.isin(missing_reacs.index)]
        # check, if any automatic gapfilling is possible
        if len(missing_reacs) == 0:
            return None
        # transform table into EC-number vs. list of NCBI protein IDs
        eccode = missing_reacs['ec-code'].apply(pd.Series).reset_index().melt(id_vars='index').dropna()[['index', 'value']].set_index('index')
        ncbiprot = missing_reacs['ncbiprotein'].apply(pd.Series).reset_index().melt(id_vars='index').dropna()[['index', 'value']].set_index('index')
        missing_reacs = pd.merge(eccode,ncbiprot,left_index=True, right_index=True).rename(columns={'value_x':'ec-code','value_y':'ncbiprotein'})
        missing_reacs = missing_reacs.groupby(missing_reacs['ec-code']).aggregate({'ncbiprotein':'unique'}).reset_index()
        
        # Step 2: map EC to reaction(s) if possible
        # -----------------------------------------
        # via MNX, BiGG, KEGG
        reacs_mapped = map_ec_to_reac(missing_reacs)
        
        # Step 3: clean and map to model reactions
        # ----------------------------------------
        # need manual curation
        self.manual_curation['reacs'] = reacs_mapped[reacs_mapped['id'].isnull()]
        # map to model reactions
        reacs_mapped = reacs_mapped[~reacs_mapped['id'].isnull()] 
        reacs_mapped['add_to_GPR'] = reacs_mapped.apply(lambda x: self._find_reac_in_model(model,x['ec-code'],x['id'],x['via']), axis=1)
                
        return reacs_mapped
    
    
# ----------------------
# Gapfilling with BioCyc
# ----------------------

class BioCycGapFiller(GapFiller):
    """
    | @TODO: Write correct doc string
    |
    |
    |
    |These three TXT files can be obtained through creating SmartTables in BioCyc 
    |and exporting these as 
    |Spreadsheets with the parameter FrameID (one & two) or Common name (three). 
    |SmartTable one: Should contain 'Accession-2' and 'Reactions of gene', 
    |SmartTable two: Should contain 'Reaction', 'Reactants of reaction', 
    |'Products of reaction', 'EC-Number', 'Reaction-Direction' and 'Spontaneous?'
    |SmartTable three: Should contain 'Compound', 'Chemical Formula' and 'InChI-Key'.

    Args:
        - GapFiller (_type_): 
            _description_
    """

    def __init__(self, model: libModel, biocyc_gene_tbl_path: str, 
                 biocyc_reacs_tbl_path: str, fasta:str) -> None:
        super().__init__()
        self._biocyc_gene_tbl = self.biocyc_gene_tbl(biocyc_gene_tbl_path) 
        # @TODO: Das ist eigtl self.full_gene_list von GapFiller!
        self._biocyc_rxn_tbl = self.biocyc_rxn_tbl(biocyc_reacs_tbl_path)
        self._model = model
        self._fasta = fasta
        self.missing_genes = None
        self.missing_reacs = None
        self.missing_metabs = None

    @property
    def biocyc_gene_tbl(self):
        return self._BioCyc_gene_tbl
    
    @biocyc_gene_tbl.setter
    # @TODO: Hier sollten wir noch diskutieren, ob Accession-2 oder Accession-1 
    # hard coded sein sollten oder nicht
    #        Dafür müssten wir nochmal ein paar entries in BioCyc dazu anschauen.
    # Locus tags in GenBank GFF file == BioCyc Accession-2 == Old locus tags in 
    # RefSeq GFF file
    # Locus tags in RefSeq GFF file == BioCyc Accession-1
    # Label in model == Locus tag from GenBank GFF file == BioCyc Accession-2
    def biocyc_gene_tbl(self, biocyc_gene_tbl_path: str) -> pd.DataFrame:
        """Parses TSV file from BioCyc to retrieve 'Accession-2' & the 
        corresponding 'Reactions of gene'
    
        Args:
           - inpath (str): 
              Path to file from BioCyc containing the 'Accession-2' to 
              'Reactions of gene' mapping
           
        Returns:
           pd.DataFrame: 
              Table containing only rows where a 'Reaction of gene' exists
        """
        # Read table
        biocyc_genes = pd.read_table(
            biocyc_gene_tbl_path, 
            usecols=['Accession-2', 'Reactions of gene'], 
            dtype=str
            )

        # Rename columns for further use
        biocyc_genes.rename(
            columns={
                'Accession-2': 'locus_tag', 'Reactions of gene': 'id'
                }, 
            inplace=True)

        # Turn empty strings into NaNs & Remove NaNs
        biocyc_genes.replace('', np.nan, inplace=True)
        biocyc_genes.dropna(inplace=True)
        
        return biocyc_genes

    @property
    def biocyc_rxn_tbl(self):
        return self._biocyc_rxn_tbl

    @biocyc_rxn_tbl.setter
    def biocyc_rxn_tbl(biocyc_reacs_tbl_path: str) -> pd.DataFrame:
        """Parses TSV file from BioCyc to retrieve 'Reaction', 'Object ID', 
        'Reactants of reaction', 'Products of reaction', 'EC-Number', 
        'Reaction-Direction' & 'Spontaneous?'

        Args:
           - biocyc_reacs_tbl_path (str):   
              Path to file from BioCyc containing the following columns:
              'Reaction' 'Object ID' 'Reactants of reaction' 
              'Products of reaction' 'EC-Number' 'Reaction-Direction' 
              'Spontaneous?'

        Returns:
           pd.DataFrame: 
              Table containing all BioCyc reactions from provided file
        """
        # Read table
        biocyc_reacs = pd.read_table(
            biocyc_reacs_tbl_path, 
            usecols=[
                'Reaction', 'Object ID', 'Reactants of reaction', 
                'Products of reaction', 'EC-Number', 'Reaction-Direction', 
                'Spontaneous?'
                ],
            dtype=str
            )

        # Rename columns for further use
        biocyc_reacs.rename(
            columns={
                'Reaction': 'equation', 'Object ID': 'id',
                'Reactants of reaction': 'Reactants', 
                'Products of reaction': 'Products', 'EC-Number': 'ec-code'
                },
            inplace=True
            )

        # Turn empty strings into NaNs
        biocyc_reacs.replace('', np.nan, inplace=True)

        # Specify empty entries in 'Spontaneous?' as False
        # @TODO Does this really make sense?
        biocyc_reacs['Spontaneous?'] = biocyc_reacs['Spontaneous?'].fillna('F')

        return biocyc_reacs

    def get_missing_genes(self):
        """Retrieves the missing genes and reactions from the BioCyc table 
        according to the 'Accession-2' identifiers
        """
        
        # Step 1: get genes from model
        # ----------------------------
        geneps_in_model = [
            _.getLabel() 
            for _ in self._model.getPlugin(0).getListOfGeneProducts()
            ]

        # Step 2: Get genes of organism from BioCyc
        # -----------------------------------------
        # For now see setter: BioCyc_gene_tbl
        
        # Step 3: BioCyc vs. model genes -> get missing genes for model
        # -------------------------------------------------------------
        self.missing_genes = self.biocyc_gene_tbl[
            ~self.biocyc_gene_tbl['locus_tag'].isin(
                geneps_in_model['locus_tag']
                )
            ]

        # Step 4: Get ncbiprotein IDs
        # ---------------------------
        # Read in FASTA file to obtain locus_tag2ncbiportein mapping
        # @TODO: GFF parsing better?
        locus_tag2ncbiprotein_df = parse_fasta_headers(self._fasta)
        locus_tag2ncbiprotein_df.rename(
            columns={'protein_id': 'ncbiprotein'},
            inplace=True
            )

        # Get the complete missing genes dataframe with the ncbiprotein IDs
        self.missing_genes = self.missing_genes.merge(
            locus_tag2ncbiprotein_df, on='locus_tag'
            )

        # Step 5: Get amount of missing genes from BioCyc for statistics
        # --------------------------------------------------------------
        self._statistics['Protein']['Total'] = len(
            self.missing_genes['locus_tag'].unique().tolist()
            )

        # Step 6: Prepare results
        # -----------------------
        # DataFrame for return
        missing_genes_res = self.missing_genes[['locus_tag', 'ncbiprotein']]
        # DataFrame for missing reactions
        self.missing_genes = self.missing_genes[['ncbiprotein', 'id']]

        return missing_genes_res

    def get_missing_reacs(self):
        """Subsets the BioCyc table with the following columns:
         'Reaction' 'Object ID' 'Reactants of reaction' 'Products of reaction' 
         'EC-Number' 'Reaction-Direction' 'Spontaneous?' to obtain the missing 
         reactions with all the corresponding data
        """
        # @TODO 
        # map to model reactions
        self._find_reac_in_model()
        
        # Step 1: get reactions from model
        # --------------------------------
        reac_model_list = get_model_reacs_or_metabs(self._model)

        # Step 2: filter missing gene list + extract ECs
        # ----------------------------------------------
        # Expand missing genes result table to merge with Biocyc reactions table
        missing_genes = pd.DataFrame(
           self.missing_genes['Reaction'].str.split('//').tolist(), 
           index=self.missing_genes['ncbiprotein']
           ).stack()
        missing_genes = missing_genes.reset_index(
            [0, 'ncbiprotein']
            )
        missing_genes.columns = ['ncbiprotein', 'Reaction']
        missing_genes['id'] = missing_genes['id'].str.strip()

        # Get missing reactions from missing genes
        missing_reacs = self.missing_genes.merge(
            self.biocyc_rxn_tbl, on='id'
            )

        # Turn entries with '//' into lists
        missing_reacs['Reactants'] = missing_reacs['Reactants'].str.split('\s*//\s*')
        missing_reacs['Products'] = missing_reacs['Products'].str.split('\s*//\s*')
        missing_reacs['ec-code'] = missing_reacs['ec-code'].str.split('\s*//\s*')

        # Step 3: Get content for column ncbiprotein
        # ------------------------------------------
        # Turn ncbiprotein column into lists of ncbiprotein IDs per reaction
        ncbiprotein_as_list = missing_reacs.groupby('id')['ncbiprotein'].apply(list).reset_index(name='ncbiprotein')
        missing_reacs.drop('ncbiprotein', axis=1, inplace=True)
        missing_reacs = ncbiprotein_as_list.merge(missing_reacs, on='id')

        # Step 4: Get amount of missing reactions from BioCyc for statistics
        # ------------------------------------------------------------------
        self._statistics['Reaction']['Total'] = len(
            missing_reacs['id'].unique().tolist()
            )

        # Step 5: Map to model reactions & cleanup
        # ----------------------------------------
        # Add column 'via'
        missing_reacs['via'] = 'BioCyc'

        # Filter reacs for already in model
        missing_reacs['add_to_GPR'] = missing_reacs.apply(
            lambda x: 
                self._find_reac_in_model(self._model,x['ec-code'],x['id'],x['via']), axis=1
            )
        
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
# # Skeleton for functions that could be used for a lab strain/organism which is in no database contained
# def get_genes_from_gff():
#     pass
# 
# 
# def get_related_metabs_reactions_blast():
#     pass
# 
# 
# def gff_gene_comp():
#     pass
# 
#

# ---------------------------------
# GapFilling with GFF and Swissprot
# ---------------------------------
class GeneGapFiller(GapFiller):
    
    GFF_COLS = {'locus_tag':'locus_tag', 
                    'eC_number':'ec-code', 
                    'protein_id':'ncbiprotein'} # :meta:
    
    def __init__(self) -> None:
        super().__init__()
        
    def get_missing_genes(self,gffpath,model:libModel) -> pd.DataFrame:
    
        # get all CDS from gff
        all_genes = parse_gff_for_cds(gffpath,self.GFF_COLS)
        # get all genes from model by locus tag
        model_locustags = [g.getLabel() for g in model.getPlugin(0).getListOfGeneProducts()]
        # filter
        missing_genes = all_genes.loc[~all_genes['locus_tag'].isin(model_locustags)]
        # formatting
        for col in self.GFF_COLS.values():
            if col not in missing_genes.columns:
                missing_genes[col] = None
                
        # save genes with no locus tag for manual curation
        self.manual_curation['gff no locus tag'] = missing_genes[missing_genes['locus_tag'].isna()]['ncbiprotein']
        
        # output
        # ncbiprotein | locus_tag | ec-code
        missing_genes =  missing_genes[~missing_genes['locus_tag'].isna()]
        missing_genes = missing_genes.explode('ncbiprotein')
        return missing_genes
    
    def get_missing_reacs(self, model:cobra.Model, 
                          missing_genes:pd.DataFrame, 
                          # prefix for pseudo ncbiprotein ids
                          prefix:str='refinegems',
                          # NCBI params
                          mail:str=None, 
                          check_NCBI:bool=False,
                          # SwissProt
                          fasta:str=None, 
                          dmnd_db:str=None, 
                          swissprot_map:str=None,
                          **kwargs) -> tuple:
        
        # Case 1:  no EC
        # --------------
        case_1 = missing_genes[missing_genes['ec-code'].isna()]
        not_case_1 = missing_genes[~missing_genes['ec-code'].isna()]
        if len(case_1) > 0:
            
            # Option 1: BLAST against SwissProt
            # +++++++++++++++++++++++++++++++++    
            # -> BLAST (DIAMOND) against SwissProt to get EC/BRENDA 
            # @TEST
            if fasta and dmnd_db and swissprot_map:
                case_1_mapped = get_ec_via_swissprot(fasta,dmnd_db,
                                            case_1,
                                            swissprot_map,
                                            **kwargs) # further optional params for the mapping
                case_1.drop('ec-code', inplace=True, axis=1)
                case_1 = case_1.merge(case_1_mapped, on='locus_tag', how='left')
                not_case_1['UniProt'] = None
                
            # Option 2: Use ML to predict EC
            # ++++++++++++++++++++++++++++++
            # @TODO
            # -> sth like DeepECTransformer (tools either not good, 
            #    not installable or no license)
            # -> use ECRECer Web service output as input 
            #    (whole protein fasta -> wait for Gwendolyn's results) -> does not work / not return
            # -> use CLEAN webservice
            #    same problem as above with the web tool

        self.manual_curation['no ncbiprotein, no EC'] = case_1[case_1['ncbiprotein'].isna() & case_1['ec-code'].isna()] 
        
        mapped_reacs = pd.concat([case_1[~(case_1['ncbiprotein'].isna() & case_1['ec-code'].isna())],not_case_1])

        # convert NaNs to None
        mapped_reacs.mask(mapped_reacs.isna(), other=None, inplace=True)

        # Case 2: still no EC but ncbiprotein
        # -----------------------------------
        #       -> access ncbi for ec (optional) 
        # @DEBUGGING ...................
        mapped_reacs = mapped_reacs.iloc[300:350,:]
        print(UserWarning('Running in debugging mode.'))
        # ..............................
        if check_NCBI and mail:
            mapped_reacs['ec-code'] = mapped_reacs.progress_apply(lambda x: get_ec_from_ncbi(mail,x['ncbiprotein']) if not x['ec-code'] and not x['ncbiprotein'].isna() else x['ec-code'], axis=1)
        
        # save entries with no EC for manual curation
        self.manual_curation['no EC'] = mapped_reacs[mapped_reacs['ec-code'].isna()]
        mapped_reacs = mapped_reacs[~mapped_reacs['ec-code'].isna()]
        
        # check, if any automatic gapfilling is still possible
        if len(mapped_reacs) == 0:
            return None
        
        # create pseudoids for entries with no ncbiprotein id
        mapped_reacs['ncbiprotein'] = mapped_reacs.apply(lambda x: f'{prefix}_{x["locus_tag"]}' if not x['ncbiprotein'] else x['ncbiprotein'], axis=1)

        # Case 3: EC found
        # ----------------
        
        # update the gene information
        updated_missing_genes = mapped_reacs.copy()
        
        # reformat missing reacs 
        mapped_reacs.drop(['UniProt','locus_tag'], inplace=True, axis=1)
        
        # transform table into EC-number vs. list of NCBI protein IDs
        # @TODO make a func out of this - occurs on multiple occasions
        eccode = mapped_reacs['ec-code'].apply(pd.Series).reset_index().melt(id_vars='index').dropna()[['index', 'value']].set_index('index')
        ncbiprot = mapped_reacs['ncbiprotein'].apply(pd.Series).reset_index().melt(id_vars='index').dropna()[['index', 'value']].set_index('index')
        mapped_reacs = pd.merge(eccode,ncbiprot,left_index=True, right_index=True).rename(columns={'value_x':'ec-code','value_y':'ncbiprotein'})
        mapped_reacs = mapped_reacs.groupby(mapped_reacs['ec-code']).aggregate({'ncbiprotein':'unique'}).reset_index()
        
        # map EC to reactions
        mapped_reacs = map_ec_to_reac(mapped_reacs[['ec-code','ncbiprotein']])
        
        # @TODO the stuff below also appear multiple times
        # save for manual curation
        self.manual_curation['reacs, no mapping'] = mapped_reacs[mapped_reacs['id'].isnull()]
        # map to model
        mapped_reacs = mapped_reacs[~mapped_reacs['id'].isnull()]
        mapped_reacs['add_to_GPR'] = mapped_reacs.apply(lambda x: self._find_reac_in_model(model,x['ec-code'],x['id'],x['via']), axis=1)
        
        return updated_missing_genes, mapped_reacs
    
    
    
    
############################################################################
# functions
############################################################################




############################################################################
# For filtering
############################################################################
# Evtl hier: 
'''
biocyc_reacs.dropna(
    subset=['id', 'Reactants', 'Products', 'EC'], inplace=True
    )
'''

''' Mapping to BiGG + addition of more information on BiGG Reactions
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
'''

# Inspired by Dr. Reihaneh Mostolizadeh's function to add BioCyc reactions to a model
def replace_reaction_direction_with_fluxes(missing_reacs: pd.DataFrame) -> pd.DataFrame:
   """Extracts the flux lower & upper bounds for each reaction through the entries in column 'Reaction-Direction'
   
   Args:
      - missing_reacs (pd.DataFrame): 
         Table containing reactions & the respective Reaction-Directions
         
   Returns:
      pd.DataFrame: 
         Input table extended with the fluxes lower & upper bounds obtained from 
         the Reaction-Directions
   """
    
   def get_fluxes(row: pd.Series) -> dict[str: str]:
      direction = row['Reaction-Direction']
      fluxes = {}
      
      if type(direction) == float:
         # Use default bounds as described in readthedocs from COBRApy
         fluxes['lower_bound'] = 'cobra_0_bound'
         fluxes['upper_bound'] = 'cobra_default_ub'
      elif 'RIGHT-TO-LEFT' in direction:
         fluxes['lower_bound'] = 'cobra_default_lb'
         fluxes['upper_bound'] = 'cobra_0_bound'
      elif 'LEFT-TO-RIGHT' in direction:
         fluxes['lower_bound'] = 'cobra_0_bound'
         fluxes['upper_bound'] = 'cobra_default_ub'
      elif 'REVERSIBLE' in direction:
         fluxes['lower_bound'] = 'cobra_default_lb'
         fluxes['upper_bound'] = 'cobra_default_ub'
      
      return str(fluxes)
   
   missing_reacs['fluxes'] = missing_reacs.apply(get_fluxes, axis=1)
   missing_reacs.drop('Reaction-Direction', axis=1, inplace=True)
   
   return missing_reacs
