# Function originally from refineGEMs.genecomp/refineGEMs.KEGG_analysis --- Modified
def compare_gene_lists(gps_in_model: pd.DataFrame, db_genes: pd.DataFrame) -> pd.DataFrame:
   
   in_db = db_genes.set_index('Locus_tag')
   in_model = gps_in_model.set_index(0)
   genes_in_db_not_in_model = in_db[~in_db.index.isin(in_model.index)]
   return genes_in_db_not_in_model.reset_index()

def get_genes_from_gff():
   pass


def compare_model_to_gff():
   pass


def get_related_metabs_reactions_blast():
   pass


def gapfill():
   ''' Main function to gapfill a model with comparison to KEGG/BioCyc/(Genbank) GFF file
   '''
   pass
