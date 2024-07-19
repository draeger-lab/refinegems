"""_summary_
"""

__author__ = "Famke Baeuerle, Gwendolyn O. Döbel and Carolin Brune"

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

from ..utility.io import parse_gff_for_gp_info
from ..curation.db_access.kegg import parse_KEGG_gene

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
        
    @abstractmethod
    def load_gene_list(self):
        pass

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
        self.report = dict()
        
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
        genes_in_model = get_model_genes(model)
        
        # Step 2: get genes of organism from KEGG
        gene_KEGG_list = KEGG().list(self.organismid)
        gene_KEGG_table = pd.read_table(io.StringIO(gene_KEGG_list), header=None)
        gene_KEGG_table.columns = ['orgid:locus','CDS','position','protein']
        gene_KEGG_table = gene_KEGG_table[['orgid:locus']]
        
        # Step 3: KEGG vs. model genes -> get missing genes for model
        genes_not_in_model = gene_KEGG_table[~gene_KEGG_table['orgid:locus'].isin(genes_in_model['orgid:locus'])]
        
        # Step 4: extract locus tag
        genes_not_in_model['locus_tag'] = genes_not_in_model['orgid:locus'].str.split(':').str[1]
        
        # Step 5: map to EC via KEGG
        geneKEGG_mapping = pd.DataFrame.from_dict(list(genes_not_in_model['orgid:locus'].apply(parse_KEGG_gene)))
        genes_not_in_model = genes_not_in_model.merge(geneKEGG_mapping, how='left', on='orgid:locus')
        
        # @TODO : What to report where and when
        self.report['missing genes (total)'] = len(genes_not_in_model)
        
        return genes_not_in_model 
    
    
    def get_missing_reacs(self, model, genes_not_in_model):
        
        # Step 1: get reactions from model 
        reac_model_list = model.getListOfReactions()
        reac_model_list = [_.id[2:] for _ in reac_model_list] # crop 'R_' prefix
        reac_model_table = pd.DataFrame({'id':reac_model_list}) # @TODO : only uses the ID
        
        # Step 2: filter missing gene list + extract ECs
        # @TODO: what should happen, if no ec-code was found -> output sth?
        reac_options = genes_not_in_model[['ec-code','ncbiprotein']]        # get relevant infos for reacs
        missing_reacs = reac_options[['ec-code','ncbiprotein']].dropna()    # drop nas
        # transform table into EC-number vs. list of NCBI protein IDs
        eccode = missing_reacs['ec-code'].apply(pd.Series).reset_index().melt(id_vars='index').dropna()[['index', 'value']].set_index('index')
        ncbiprot = missing_reacs['ncbiprotein'].apply(pd.Series).reset_index().melt(id_vars='index').dropna()[['index', 'value']].set_index('index')
        missing_reacs = pd.merge(eccode,ncbiprot,left_index=True, right_index=True).rename(columns={'value_x':'ec-code','value_y':'ncbiprotein'})
        missing_reacs.groupby(missing_reacs['ec-code']).aggregate({'ncbiprotein':'unique'}).reset_index()
                
        # Step 3: mapping based on EC number (via KEGG)
        
        # Step 4: map to BiGG

        # Step 6: compare to model 
        
        pass 
    
    
# ----------------------
# Gapfilling with BioCyc
# ----------------------

class BioCycGapFiller(GapFiller):

    def __init__(self, model: libModel, biocyc_gene_tbl_path: str, biocyc_reacs_tbl_path: str) -> None:
        super().__init__()
        self._biocyc_gene_tbl = self.biocyc_gene_tbl(biocyc_gene_tbl_path) # @TODO: Das ist eigtl self.full_gene_list von GapFiller!
        self._biocyc_reacs_tbl = self.biocyc_reacs_tbl(biocyc_reacs_tbl_path)
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
    def biocyc_reacs_tbl(self):
        return self._biocyc_reacs_tbl

    @biocyc_reacs_tbl.setter
    def biocyc_reacs_tbl(biocyc_reacs_tbl_path: str) -> pd.DataFrame:
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
                            {'Reactants of reaction': 'Reactants', 'Products of reaction': 'Products', 'EC-Number': 'EC'},
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
        geneps_in_model = [_.getLabel() for _ in self.model.getPlugin(0).getListOfGeneProducts()]

        # Step 2: Get genes of organism from BioCyc
        # For now see setter: BioCyc_gene_tbl
        
        # Step 3: BioCyc vs. model genes -> get missing genes for model
        self.missing_biocyc_genes = self.biocyc_gene_tbl[~self.biocyc_gene_tbl['locus_tag'].isin(geneps_in_model['locus_tag'])]

        # Step 4: Get amount of missing genes from BioCyc for statistics
        self._statistics['Protein']['Total'] = len(self.missing_biocyc_genes['locus_tag'].unique().tolist())

    def get_missing_reacs(self) -> tuple[tuple[pd.DataFrame, pd.DataFrame], pd.DataFrame]:
        """Subsets the BioCyc table with the following columns: 
         'Reaction' 'Reactants of reaction' 'Products of reaction' 'EC-Number' 'Reaction-Direction' 'Spontaneous?'
         to obtain the missing reactions with all the corresponding data 
         & Adds the according BiGG Reaction identifiers

        Args:
           - inpath (str): 
                Path to file from BioCyc containing the following columns:
                'Reaction' 'Reactants of reaction' 'Products of reaction' 'EC-Number' 
                'Reaction-Direction' 'Spontaneous?'
      
        Returns:
           tuple: 
                Tuple (1) & table (2)

                (1) tuple: 
                        Two tables (1) & (2)

                        (1) pd.DataFrame: Table containing only the metabolites corresponding to the missing BioCyc reactions
                        (2) pd.DataFrame: Table containing only the metabolites corresponding to the missing BiGG reactions
         
                (2) pd.DataFrame: 
                        Table containing the missing reactions with the corresponding data
        """
   
        missing_biocyc_reactions = pd.DataFrame(
           self.missing_biocyc_genes['Reaction'].str.split('//').tolist(), index=self.missing_biocyc_genes['locus_tag']
           ).stack()
        missing_biocyc_reactions = missing_biocyc_reactions.reset_index([0, 'locus_tag'])
        missing_biocyc_reactions.columns = ['locus_tag', 'Reaction']
        missing_biocyc_reactions['Reaction'] = missing_biocyc_reactions['Reaction'].str.strip()

        model_reacs = get_model_reacs_or_metabs(model_libsbml)
        biocyc_reacs = get_biocyc_reactions(inpath)

        # Get missing reactions from missing genes
        missing_reactions = genes2reaction.merge(biocyc_reacs, on='Reaction')

        # Turn entries with '//' into lists
        missing_reactions['Reactants'] = missing_reactions['Reactants'].str.split('\s*//\s*')
        missing_reactions['Products'] = missing_reactions['Products'].str.split('\s*//\s*')
        missing_reactions['EC'] = missing_reactions['EC'].str.split('\s*//\s*')

        # Turn locus_tag column into lists of locus tags per reaction
        locus_tags_as_list = missing_reactions.groupby('Reaction')['locus_tag'].apply(list).reset_index(name='locus_tag')
        missing_reactions.drop('locus_tag', axis=1, inplace=True)
        missing_reactions = locus_tags_as_list.merge(missing_reactions, on='Reaction')
        self._statistics['Reaction']['Total'] = len(missing_reactions['Reaction'].unique().tolist())

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
        return (biocyc_metabs_from_reacs, bigg_metabs_from_reacs), missing_reactions 

   

   
        return missing_biocyc_reactions

    
# ----------------
# Gapfilling no DB
# ----------------

class GeneGapFiller(GapFiller):
    pass

############################################################################
# functions
############################################################################