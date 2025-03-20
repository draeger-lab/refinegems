import pandas as pd

from libsbml import ListOfGeneProducts


def extend_gp_annots_via_files(gene_list: ListOfGeneProducts, ids2annots: str, email: str, gff_paths: list[str], protein_fasta: str, filename: str, lab_strain: bool=False) -> None:

    # 1. If no table available

        # 1. (optional) Get information from GFF, if both from both, if none then ignore, if one then from that one
        # -> Table with locus_tag, Protein ID, RefSeq

        # 2. If Protein ID/RefSeq/both is empty extract ID from model ID (Extend table with model_id; also add column id_from_model_id?)

        # 3. Extend table with name -> Query NCBI for name via
        #   1. Protein ID, if available & RefSeq ID failed
        #   2. RefSeq ID, if available & Protein ID failed
        #   3. model_id, if both of the above available & failed

    # 2. Use table to fill in information in model
    
    pass