import pandas as pd

from cobra.io.sbml import _f_gene
from functools import reduce
from libsbml import ListOfGeneProducts

from refinegems.utility.io import get_gff_variety, parse_gff_for_cds


def get_gpid2name_mapping(gene_list: ListOfGeneProducts, gff_paths: list[str], email: str) -> pd.DataFrame:

    # 1. Get model IDs & potential contained IDs
    print('Extracting model IDs and potential valid database IDs from model...')
    modelid2potetialid = {'model_id': [], 'database_id': []}
    for gene in gene_list:
        # Get current gene id
        gene_id = gene.getId()
        modelid2potetialid['model_id'].append(gene_id)

        # Get gene_id without prefix & replace ascii numbers with according characters
        gene_id = _f_gene(gene_id)

        # Extract potential database ID from gene_id
        if (gene_id[0] == 'W'): #addition to work with KC-Na-01
            # @TODO Implement better way of handling this -> Handle with COBRApy
            entry = gene_id[:-7] if '__46__' in gene_id else gene_id[:-2] # Required for VMH models
            name, locus = search_ncbi_for_gpr(entry)
            gene.setName(name)
            if (locus2id is not None) and (entry in locus2id.index):
                locus = locus2id.loc[entry, 'LocusTag']
        
        elif (gene_id != 'spontaneous') and (gene_id != 'Unknown'): # Has to be omitted as no additional data can be retrieved neither from NCBI nor the CarveMe input file
            if 'prot_' in gene_id:
                id_string = gene_id.split(r'prot_')[1].split(r'_')  # All NCBI CDS protein FASTA files have the NCBI protein identifier after 'prot_' in the FASTA identifier
                ncbi_id = id_string[0]  # If identifier contains no '_', this is full identifier
            else:
                id_string = gene_id.split(r'_')
                if 'peg' in id_string: 
                    pass
              
            if len(id_string) == 2: # Can be the case if ID is locus tag, for example
                continue # Ignore locus tags as no valid identifiers
            if (len(id_string) > 2):  # Identifier contains '_'
            # Check that the second entry consists of a sequence of numbers -> Valid RefSeq identifier! 
            # (Needs to be changed if there are other gene idenitfiers used that could contain '_' & need to be handled differently)
                if re.fullmatch(r'^\d+\d+$', id_string[1], re.IGNORECASE):
                    # Merge the first two parts with '_' as this is complete identifier
                    # Merge the resulting string with the third string in i_string to get complete identifier with version spec
                    ncbi_id = f'{"_".join(id_string[:2])}.{id_string[2]}'

            # If identifier matches RefSeq ID pattern   
            if re.fullmatch(r'^(((AC|AP|NC|NG|NM|NP|NR|NT|NW|WP|XM|XP|XR|YP|ZP)_\d+)|(NZ_[A-Z]{2,4}\d+))(\.\d+)?$', ncbi_id, re.IGNORECASE):
                name, locus = search_ncbi_for_gpr(ncbi_id)
                if (locus2id is not None) and (ncbi_id in locus2id.index):
                    locus = locus2id.loc[ncbi_id, 'LocusTag']

            # If identifier only contains numbers 
            # -> Get the corresponding data from the CarveMe input file
            elif re.fullmatch(r'^\d+$', ncbi_id, re.IGNORECASE):
                if id2locus_name is not None:
                    name, locus = id2locus_name[id2locus_name['protein_id']==ncbi_id][['name', 'locus_tag']].values[0]
                else: 
                    genes_missing_annotation.append(gene_id)
        
            # If identifier matches ncbiprotein ID pattern
            elif re.fullmatch(r'^(\w+\d+(\.\d+)?)|(NP_\d+)$', ncbi_id, re.IGNORECASE):
                name, locus = search_ncbi_for_gpr(ncbi_id)
        
            # For lab strains use the locus tag from the annotation file   
            if lab_strain and id2locus_name is not None:
                locus = id2locus_name[id2locus_name['protein_id']==ncbi_id][['locus_tag']].values[0][0]
        
            if ncbi_id not in genes_missing_annotation:
                if name and locus:
                    gene.setName(name)
                    if len(locus) == 1: gene.setLabel(locus[0])
        # Case 1: CarveMe version; Contains NCBI Protein ID, lcl_CP035291_1_prot_QCY37541_1_361

        # Case 2: RefSeq contained without any surrounding stuff, WP_011275285.1, From CarveMe: YP_005228568_1

        # Case 3: Rast contained without any surrounding stuff, 1134914.3.1020_56.peg

        # Case 4: KEGG contained without any surrounding stuff, sha:SH1836
        
        # Case 5: Model ID contains locus_tag, AB-1_S128_00983


    # 2.A (Optional) Get information from GFF(s)
    if gff_paths:
        # Get attributes to keep per GFF variety
        gff_attributes = {
            'refseq': {'locus_tag': 'locus_tag', 'protein_id': 'REFSEQ', 'product': 'name'}, 
            'genbank': {'locus_tag': 'locus_tag', 'protein_id': 'NCBI', 'product': 'name'}, 
            'prokka': {'locus_tag': 'locus_tag', 'product': 'name'}
            }

        # Parse & collect all provided GFFs
        gffs = []
        for gffp in gff_paths:
            current_gff_variety = get_gff_variety(gffp)
            current_gff = parse_gff_for_cds(gffp, keep_attributes=gff_attributes[current_gff_variety])
            gffs.append(current_gff)

        # Merge all GFFs on common column locus_tag
        if len(gffs) == 1:
            gff_mapping = gffs[0]
        else:
            gff_mapping = reduce(
                lambda left, right: pd.merge(left, right, on=['locus_tag'], how='outer', sort=True), gffs
                )


def extend_gp_annots_via_files(gene_list: ListOfGeneProducts, gff_paths: list[str], email: str, filename: str, lab_strain: bool=False) -> None:

    # 1. If no table available



        # 2. If Protein ID/RefSeq/both is empty extract ID from model ID (Extend table with model_id; also add column id_from_model_id?)

        # 3. Extend table with name -> Query NCBI for name via
        #   1. Protein ID, if available & RefSeq ID failed
        #   2. RefSeq ID, if available & Protein ID failed
        #   3. model_id, if both of the above available & failed

    # 2. Use table to fill in information in model
    
    pass