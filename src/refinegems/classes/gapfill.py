#!/usr/bin/env python
"""Add reactions, genes and more to a model based on different gap-filling methods.
All (current) algorithms are separated into three steps: finding missing genes, finding missing reactions
and trying to add the found as missing entities to the model.

Available gap filling methods:

- | :py:class:`~refinegems.classes.gapfill.KEGGapFiller`
  | Mainly utilises information from the KEGG database. Needs a KEGG organism ID.
  | Estimated runtime: *to be determined*

- | :py:class:`~refinegems.classes.gapfill.BioCycGapFiller`
  | Mainly utilises information from the BioCyc database. Requires access to BioCyc SmartTables.
  | Estimated runtime: *to be determined*

- | :py:class:`~refinegems.classes.gapfill.GeneGapFiller`
  | Search for gaps using the GFF file and information from SwissProt.
  | Estimated runtime: *to be determined*

"""

__author__ = "Famke Baeuerle, Gwendolyn O. Döbel, Carolin Brune, Finn Mier and Dr. Reihaneh Mostolizadeh"

############################################################################
# requirements
############################################################################

from abc import ABC, abstractmethod

import cobra
from cobra.io.sbml import _f_gene, _f_gene_rev
from copy import deepcopy
import io
import json
import logging
import math
import numpy as np
import os
import pandas as pd
import re
import warnings

from bioservices.kegg import KEGG
from itertools import chain
from libsbml import Model as libModel
from pathlib import Path
from tempfile import NamedTemporaryFile
from typing import Literal, Union

from tqdm import tqdm

tqdm.pandas()

from ..developement.decorators import *
from ..utility.db_access import (
    get_ec_from_ncbi,
    get_ec_via_swissprot,
    parse_KEGG_gene,
    parse_KEGG_ec,
    map_to_homologs,
)
from ..utility.io import (
    load_a_table_from_database,
    parse_gff_for_cds,
    load_model,
    write_model_to_file,
)
from ..utility.entities import (
    create_gp,
    create_gpr,
    build_reaction_bigg,
    build_reaction_kegg,
    build_reaction_mnx,
    isreaction_complete,
    parse_reac_str,
    validate_reaction_compartment_bigg,
    REF_COL_GF_GENE_MAP,
)
from .medium import Medium, medium_to_model
from .reports import GapFillerReport

# @Note:
#   some reactions have @DEBUGGING
#   enable these parts to shorten runtime during debugging (used to work on a subset,
#   not on the whole input)

############################################################################
# variables
############################################################################


############################################################################
# functions
############################################################################


# Mapping for BioCyc Reactions
# ----------------------------
def map_biocyc_to_reac(
    biocyc_reacs: pd.DataFrame, use_MNX: bool = True, use_BiGG: bool = True
) -> pd.DataFrame:
    """Based on a table containing BioCyc reactions,
    map them to reactions in other databases (if a mapping is possible)

    Args:
        - biocyc_reacs (pd.DataFrame):
            Table containing BioCyc reactions information.
            Should contain the columns: Reaction | Object ID | EC-Number | Spontaneous?
            Can be doenloaded as a SmartTable from BioCyc.
        - use_MNX (bool, optional):
            Try mapping using the MetaNetX database.
            Defaults to True.
        - use_BiGG (bool, optional):
            Try mapping using the BiGG database.
            Defaults to True.

    Returns:
        pd.DataFrame:
            The extended table.
    """

    def _clean_res_row(
        res_row: pd.Series, mapped2db: Literal["BiGG", "MetaNetX"]
    ) -> pd.Series:
        """Clean a row of a table containing mapping results from BioCyc to another database

        Args:
            - res_row (pd.Series):
                Row containing a mapping from BioCyc to another database
            - mapped2db (Literal['BiGG', 'MetaNetX']):
                Database where Biocyc entries where mapped to. One of 'BiGG', 'MetaNetX'.

        Returns:
            pd.Series:
                The cleaned row.
        """
        # Mapping for equation column name for other databases
        dbeq2eq = {"BiGG": "reaction_string", "MetaNetX": "mnx_equation"}

        # Get list of mapped to database IDs in row if available
        id_list = list(res_row[mapped2db]) if res_row[mapped2db] else None
        if id_list:
            # Move BioCyc IDs to references
            biocyc_reac_id = res_row["id"]
            res_row["id"] = str(id_list[0])  # Ensure ID is a string

            # Remove ID from column 'id' from list in column 'alias'
            if len(id_list) != 1:
                id_list.remove(res_row["id"])
                alias = id_list
            else:
                alias = None

            # Move BioCyc IDs and alias to references
            res_row["reference"] = {"metacyc.reaction": biocyc_reac_id, "alias": alias}

            # Move equation from other database to equation column
            res_row["equation"] = res_row[dbeq2eq.get(mapped2db)]

            # Replace BioCyc in via column with mapped to database
            res_row["via"] = mapped2db

        return res_row

    def _map_biocyc_to_mnx(unmapped_reacs: pd.DataFrame) -> pd.DataFrame:
        """
        | Helper function for :py:func:`~refinegems.classes.gapfill.map_biocyc_to_reac`
        | Maps Biocyc IDs in a table to the corresponding MetaNetX IDs &
        | Adds information obtained via the MetaNetX IDs to the resulting table

        Args:
            - unmapped_reacs (pd.DataFrame):
                Table containing BioCyc IDs

        Returns:
            pd.DataFrame:
                The extended table.
        """

        # Step 1: Mapping & Obtain information
        # ------------------------------------
        # Mapping to MetaNetX + addition of more information on MetaNetX
        # Reactions
        # Get MetaNetX BioCyc
        mnx2biocyc_reacs = load_a_table_from_database(
            """
            SELECT x.source, x.id, p.mnx_equation, p.\'ec-code\' 
            FROM mnx_reac_xref x INNER JOIN mnx_reac_prop p 
            USING(id) 
            WHERE x.source LIKE \'%metacyc%\' 
            OR x.source LIKE \'%biocyc%\'
            """
        )

        # Change column names to fit input table
        mnx2biocyc_reacs.rename(
            {"id": "MetaNetX", "source": "id"}, axis=1, inplace=True
        )

        # Transform id column to only contain BioCyc/MetaCyc IDs
        mnx2biocyc_reacs["id"] = mnx2biocyc_reacs["id"].str.split(":", n=1).str[-1]

        # Crop table to contain MetaNetX IDs set per BioCyc ID
        mnx_as_list = (
            mnx2biocyc_reacs.groupby("id")["MetaNetX"]
            .apply(set)
            .reset_index(name="MetaNetX")
        )
        mnx2biocyc_reacs.drop("MetaNetX", axis=1, inplace=True)
        mnx2biocyc_reacs = mnx_as_list.merge(mnx2biocyc_reacs, on="id")

        # Drop duplicates to get the unique BioCyc IDs
        mnx2biocyc_reacs.drop_duplicates(subset="id", inplace=True)

        # Merge mnx2biocyc_reacs with unmapped_reacs to get new information
        reacs_mapped = mnx2biocyc_reacs.merge(unmapped_reacs, on="id", how="right")

        # Step 2: Clean-up
        # ----------------
        # Remove all NaNs from DataFrame
        reacs_mapped.replace(np.nan, None, inplace=True)

        # Turn MetaNetX column in single value, rename column to id
        # & if multiple MetaNetX IDs exist add to column alias
        reacs_mapped = reacs_mapped.apply(_clean_res_row, args=("MetaNetX",), axis=1)

        # Create list of EC codes in column ec-code_x,
        # Join both ec-code columns into one & Create a set of ec-codes
        reacs_mapped["ec-code_x"] = reacs_mapped["ec-code_x"].str.split(r"\s*;\s*")
        reacs_mapped["ec-code_x"] = (
            reacs_mapped["ec-code_x"]
            .fillna({i: [] for i in reacs_mapped.index})
            .map(set)
            .map(list)
        )
        reacs_mapped["ec-code"] = (
            (reacs_mapped["ec-code_x"] + reacs_mapped["ec-code_y"]).map(set).map(list)
        )

        # Drop all unnecessary columns
        reacs_mapped.drop(
            ["MetaNetX", "mnx_equation", "ec-code_x", "ec-code_y"], axis=1, inplace=True
        )

        # Step 3: Return result
        # ---------------------
        return reacs_mapped

    def _map_biocyc_to_bigg(unmapped_reacs: pd.DataFrame) -> pd.DataFrame:
        """
        | Helper function for :py:func:`~refinegems.classes.gapfill.map_biocyc_to_reac`
        | Maps Biocyc IDs in a table to the corresponding BiGG IDs &
        | Adds information obtained via the BiGG IDs to the resulting table

        Args:
            - unmapped_reacs (pd.DataFrame):
                Table containing BioCyc IDs

        Returns:
            pd.DataFrame:
                The extended table.
        """

        # Step 1: Mapping & Obtain information
        # ------------------------------------
        # Mapping to BiGG + addition of more information on BiGG Reactions
        # Get BiGG BioCyc
        bigg2biocyc_reacs = load_a_table_from_database(
            "SELECT id, BioCyc, name, reaction_string from bigg_reactions"
        )

        # Filter BiGG table to contain only reaction IDs with valid compartments
        mask = bigg2biocyc_reacs.apply(
            lambda x: validate_reaction_compartment_bigg(
                parse_reac_str(x["reaction_string"], "BiGG")[2]
            ),
            axis=1,
        )
        mask.replace("exchange", True, inplace=True)
        bigg2biocyc_reacs = bigg2biocyc_reacs[mask]

        # Change column names to fit input table
        bigg2biocyc_reacs.rename({"id": "BiGG", "BioCyc": "id"}, axis=1, inplace=True)

        # Crop table to contain BiGG IDs set per BioCyc ID
        bigg_as_list = (
            bigg2biocyc_reacs.groupby("id")["BiGG"].apply(set).reset_index(name="BiGG")
        )
        bigg2biocyc_reacs.drop("BiGG", axis=1, inplace=True)
        bigg2biocyc_reacs = bigg_as_list.merge(bigg2biocyc_reacs, on="id")

        # Drop duplicates to get the unique BioCyc IDs
        bigg2biocyc_reacs.drop_duplicates(subset="id", inplace=True)

        # Merge bigg2biocyc_reacs with unmapped_reacs to get new information
        reacs_mapped = bigg2biocyc_reacs.merge(unmapped_reacs, on="id", how="right")

        # Step 2: Clean-up
        # ----------------
        # Remove all NaNs from DataFrame
        reacs_mapped.replace(np.nan, None, inplace=True)

        # Turn BiGG column in single value, rename column to id
        # & if multiple BiGG IDs exist add to column alias
        reacs_mapped = reacs_mapped.apply(_clean_res_row, args=("BiGG",), axis=1)

        # Drop all unnecessary columns
        reacs_mapped.drop(["BiGG", "reaction_string", "name"], axis=1, inplace=True)

        # Step 3: Return result
        # ---------------------
        return reacs_mapped

    # Step 1: Get missing reactions without GPR
    # -----------------------------------------
    mask = biocyc_reacs["add_to_GPR"].isna()
    biocyc_reacs_to_map = biocyc_reacs[mask]

    # Step 2: Mapping
    # ---------------
    # Map to MetaNetX
    if use_MNX:
        biocyc_reacs_to_map = _map_biocyc_to_mnx(biocyc_reacs_to_map)

    # Map to BiGG
    if use_BiGG:
        to_map = biocyc_reacs_to_map[biocyc_reacs_to_map["via"] == "BioCyc"]
        if len(to_map) > 0:
            to_map = _map_biocyc_to_bigg(to_map)
            biocyc_reacs_to_map = pd.concat(
                [to_map, biocyc_reacs_to_map[~(biocyc_reacs_to_map["via"] == "BioCyc")]]
            )

    # Step 3: Gather result(s)
    # ------------------------
    # Merge all tables if necessary
    biocyc_reacs = pd.concat(
        [biocyc_reacs_to_map, biocyc_reacs[~biocyc_reacs["add_to_GPR"].isna()]],
        sort=True,
        ignore_index=True,
    )

    # Return results
    return biocyc_reacs


# Mapping of EC numbers
# ---------------------


def map_ec_to_reac(
    table: pd.DataFrame,
    use_MNX: bool = True,
    use_BiGG: bool = True,
    use_KEGG: bool = True,
    threshold_add_reacs: int = 5,
) -> pd.DataFrame:
    """Based on a table of NCBI protein IDs and EC numbers,
    map them to reactions via different databases
    (if a mapping is possible).

    input table should have format:
        ``ec-code | ncbiprotein``

    output table has the format:
        ``ec-code | ncbiprotein | id | equation | reference | is_transport | via``

    Args:
        - table (pd.DataFrame):
            The input table.
        - use_MNX (bool, optional):
            Try mapping using the MetaNetX database.
            Defaults to True.
        - use_BiGG (bool, optional):
            Try mapping using the BiGG database.
            Defaults to True.
        - use_KEGG (bool, optional):
            Try mapping using the KEGG database.
            Defaults to True.
        - threshold_add_reacs (int, optional):
            Maximum number of reactions allowed per EC number per ncbiprotein ID
            to be added. Otherwise skip addition of reactions due to insufficient evidence
            Defaults to 5.

    Returns:
        pd.DataFrame:
            The extended table.
    """

    def _map_ec_to_reac_mnx(
        unmapped_reacs: pd.DataFrame, threshold_add_reacs: int = 5
    ) -> pd.DataFrame:
        """Helper function of :py:func:`~refinegems.classes.gapfill.map_ec_to_reac`
        for the mapping using the MetaNetX database.

        Args:
            - unmapped_reacs (pd.DataFrame):
                The input table (ec-code and ncbiprotein columns)
            - threshold_add_reacs (int, optional):
                Maximum number of reactions allowed per EC number per ncbiprotein ID
                to be added. Otherwise skip addition of reactions due to insufficient evidence
                Defaults to 5.

        Returns:
            pd.DataFrame:
                The extended table
        """

        # input: pd.DataFrame with at least the ec-code column
        # load MNX reac prop table
        mnx_reac_prop = load_a_table_from_database("mnx_reac_prop", False)
        # convert table into one EC-number per row
        mnx_reac_prop.drop("is_balanced", inplace=True, axis=1)
        mnx_reac_prop["ec-code"] = mnx_reac_prop["ec-code"].apply(
            lambda x: x.split(r";") if isinstance(x, str) else None
        )
        # exclude entries without EC-number
        mnx_reac_prop = mnx_reac_prop.explode("ec-code").dropna(subset="ec-code")
        # merge with unmapped reactions
        reacs_mapped = unmapped_reacs.merge(mnx_reac_prop, on="ec-code", how="left")
        # filter out mappings with more than x merged reactions, to ensure quality
        reacs_mapped = reacs_mapped.explode("ncbiprotein")
        reacs_mapped = reacs_mapped.groupby(["ncbiprotein", "ec-code"]).filter(
            lambda x: len(x) < threshold_add_reacs
        )
        reacs_mapped = reacs_mapped.merge(
            unmapped_reacs.explode("ncbiprotein"),
            on=["ec-code", "ncbiprotein"],
            how="outer",
        )
        # rename columns and cleanup
        reacs_mapped.rename({"mnx_equation": "equation"}, inplace=True, axis=1)
        reacs_mapped = (
            reacs_mapped.groupby(["ec-code", "id"])
            .agg(
                {
                    "ncbiprotein": lambda x: x.tolist(),
                    "equation": "first",
                    "reference": "first",
                    "is_transport": "first",
                }
            )
            .reset_index()
        )
        reacs_mapped = reacs_mapped.reindex(
            columns=[
                "ec-code",
                "ncbiprotein",
                "id",
                "equation",
                "reference",
                "is_transport",
            ]
        )
        reacs_mapped["via"] = reacs_mapped["id"].apply(
            lambda x: "MetaNetX" if x else None
        )

        return reacs_mapped

    def _map_ec_to_reac_bigg(unmapped_reacs: pd.DataFrame) -> pd.DataFrame:
        """Helper function of :py:func:`~refinegems.classes.gapfill.map_ec_to_reac`
        for the mapping using the BiGG database.

        Args:
            - unmapped_reacs (pd.DataFrame):
                The input table (ec-code and ncbiprotein columns)

        Returns:
            pd.DataFrame:
                The extended table
        """

        # load BiGG reaction namespace
        bigg_reacs = load_a_table_from_database("bigg_reactions", False)
        bigg_reacs.dropna(subset="EC Number", inplace=True)
        bigg_reacs = bigg_reacs[["id", "reaction_string", "EC Number"]].rename(
            {"reaction_string": "equation", "EC Number": "ec-code"},
            inplace=False,
            axis=1,
        )
        bigg_reacs["ec-code"] = bigg_reacs["ec-code"].apply(
            lambda x: x.split(r",") if isinstance(x, str) else None
        )
        bigg_reacs = bigg_reacs.explode("ec-code")
        # merge with unmapped reactions
        bigg_mapping = unmapped_reacs.merge(bigg_reacs, on=["ec-code"], how="left")
        bigg_mapping.mask(bigg_mapping.isna(), other=None, inplace=True)
        # make conform to format
        bigg_mapping["reference"] = None
        bigg_mapping["is_transport"] = None
        bigg_mapping["via"] = bigg_mapping["id"].apply(lambda x: "BiGG" if x else None)

        return bigg_mapping

    def _map_ec_to_reac_kegg(unmapped_reacs: pd.DataFrame) -> pd.DataFrame:
        """Helper function of :py:func:`~refinegems.classes.gapfill.map_ec_to_reac`
        for the mapping using the KEGG database.

        Args:
            - unmapped_reacs (pd.DataFrame):
                The input table (ec-code and ncbiprotein columns)

        Returns:
            pd.DataFrame:
                The extended table
        """

        # get KEGG EC number information
        eccodes = unmapped_reacs["ec-code"]
        kegg_mapped = pd.DataFrame.from_dict(
            list(eccodes.progress_apply(parse_KEGG_ec))
        )
        kegg_mapped = unmapped_reacs.merge(kegg_mapped, on="ec-code")

        kegg_mapped["is_transport"] = None
        kegg_mapped["via"] = kegg_mapped["id"].apply(lambda x: "KEGG" if x else None)
        kegg_mapped.explode(column="id")

        return kegg_mapped

    # input table should have format:
    #   ec-code | ncbiprotein
    #   one EC number per row, list of ncbiprotein per row allowed
    if (
        len(table.columns) != 2
        or "ec-code" not in table.columns
        or "ncbiprotein" not in table.columns
    ):
        raise ValueError("Wrong table format. Cannot map EC to reaction.")

    # map to MetaNetX
    if use_MNX:
        table = _map_ec_to_reac_mnx(table, threshold_add_reacs)
    # map to BiGG
    if use_BiGG:
        if "id" in table.columns:
            to_map = table[table["id"].isna()][["ec-code", "ncbiprotein"]]
            if len(to_map) > 0:
                to_map = _map_ec_to_reac_bigg(to_map)
                table = pd.concat([to_map, table[~table["id"].isna()]])
        else:
            table = _map_ec_to_reac_bigg(table)

    # map to KEGG
    if use_KEGG:
        if "id" in table.columns:
            to_map = table[table["id"].isna()][["ec-code", "ncbiprotein"]]
            if len(to_map) > 0:
                to_map = _map_ec_to_reac_kegg(to_map)
                table = pd.concat([to_map, table[~table["id"].isna()]])
        else:
            table = _map_ec_to_reac_kegg(table)
    # explode
    table = table.explode("id", ignore_index=True).explode("id", ignore_index=True)

    # output:   ec-code ncbiprotein	id	equation	reference	is_transport	via
    return table


############################################################################
# classes
############################################################################

# -------------
# abtract class
# -------------


class GapFiller(ABC):
    """Abstract base class for the gap filling.

    Already includes functions for the "filling" part of the gap-filling approach
    and some helper functions. Each subclass needs an implementation of `find_missing_genes`
    and `find_missing_reactions` to determine the entities, that are missing in the model.

    Attributes:
        - full_gene_list (list):
            List of all the genes.
        - missing_genes:
            DataFrame containing all genes found as missing with additional information. <br>
            ``ncbiprotein | locus_tag | <optional columns>``
        - missing_reactions:
            DataFrame containing all reactions found as missing with additional information.<br>
            ``ec-code | ncbiprotein | id | equation | reference | <is_transport> | via | add_to_GPR``
        - geneid_type (str):
            What type of gene ID the model contains.
            Defaults to 'ncbi'.
        - _statistics (dict):
            Dictionary of statistical information of the gap-filling run. Includes e.g.
            the number of added genes and reactions.
        - manual_curation (dict):
            Dictionary of reaction and gene IDs to be used for manual curation.
    """

    def __init__(self) -> None:

        # data
        self.full_gene_list = None
        self.missing_genes = (
            None  # missing genes, that have not yet been sorted into any category
        )
        self.missing_reactions = (
            None  # missing reacs, that have not yet been sorted into any category
        )

        # general information
        self.geneid_type = "ncbi"
        self._variety = "Undefined"  # Specifies the variety of the gapfiller, e.g. 'BioCyc', 'KEGG', 'GFF + SwissProt'

        # collect stats & Co, can be extended by subclasses
        self._statistics = {
            "genes": {
                "missing (total)": 0,
                "unmappable": 0,
                "missing (mappable)": 0,
                "duplicated": 0,
                "added": 0,
                "building failed": 0,
                "missing (remaining)": 0,
            },
            "reactions": {
                "missing (based on genes)": 0,
                "already in model": 0,
                "missing (total)": 0,
                "mapped to MNX": 0,
                "mapped to BiGG": 0,
                "mapped to KEGG": 0,
                "unmappable": 0,
                "added": 0,
                "building failed": 0,
                "missing (remaining)": 0,
            },
        }
        self.manual_curation = {"genes": {}, "reactions": {}}

    # abstract methods
    # ----------------

    @abstractmethod
    def find_missing_genes(self, model):
        """Find missing genes in the model. Parameters can be extended as needed.

        Needs to save a table in the format

        ``ncbiprotein | locus_tag | <optional columns>``

        to the attribute missing_genes.

        """
        pass

    @abstractmethod
    def find_missing_reactions(self, model):
        """Find missing reactions in the model. Parameters can be extended as needed.

        Needs to save a table of the format

        ``ec-code | ncbiprotein | id | equation | reference | <is_transport> | via | add_to_GPR``

        to the attribute missing_reactions.

        Method specific information can be added to the reference column, which is
        expected to contain a dictionary. The 'via' column describes the way
        the database will be added to the model.
        """
        pass

    # finding the gaps
    # ----------------
    def _find_reac_in_model(
        self,
        model: cobra.Model,
        eccode: str,
        id: str,
        idtype: Literal["MetaNetX", "KEGG", "BiGG", "BioCyc"],
        include_ec_match: bool = False,
    ) -> Union[None, list]:
        """Helper function of :py:class:`~refinegems.classes.gapfill.GapFiller`.
        Search the model for an ID (and optionally EC number), to determine, if the
        reaction is in the model.

        Args:
            - model (cobra.Model):
                The model loaded with COBRApy.
            - eccode (str):
                The EC number in the format: X.X.X.X
            - id (str):
                The ID to search for.
            - idtype (Literal['MetaNetX','KEGG','BiGG', 'BioCyc']):
                Name of the database the ID belongs to.
            - include_ec_match (bool, optional):
                Option to include a match if only the EC number matches.
                Defaults to False.

        Returns:
            (1) Case one or more match found:
                    list:
                        List of the ID of the reactions in the model, that match
                        the query.

            (2) Case no match:
                    None:
                        Nothing found.
        """

        MAPPING = {
            "MetaNetX": "metanetx.reaction",
            "KEGG": "kegg.reaction",
            "BiGG": "bigg.reaction",
            "BioCyc": "metacyc.reaction",
        }

        found = []
        for r in model.reactions:
            if MAPPING[idtype] in r.annotation.keys():
                if (
                    isinstance(r.annotation[MAPPING[idtype]], list)
                    and id in r.annotation[MAPPING[idtype]]
                ):
                    found.append(r.id)
                elif (
                    isinstance(r.annotation[MAPPING[idtype]], str)
                    and id == r.annotation[MAPPING[idtype]]
                ):
                    found.append(r.id)
            if include_ec_match and eccode and "ec-code" in r.annotation.keys():
                if (
                    isinstance(r.annotation["ec-code"], list)
                    and eccode in r.annotation["ec-code"]
                ):
                    found.append(r.id)
                elif (
                    isinstance(r.annotation["ec-code"], str)
                    and eccode == r.annotation["ec-code"]
                ):
                    found.append(r.id)

        found = list(set(found))

        if len(found) > 0:
            return found

        return None

    # actual "Filling" part
    # ---------------------
    def add_genes_from_table(self, model: libModel, gene_table: pd.DataFrame) -> None:
        """Create new GeneProduct for a table of genes in the format:

        | ncbiprotein | locus_tag | UniProt | ... |

        The dots symbolise additional columns, that can be passed to the function,
        but will not be used by it. The other columns, except UniProt, are required.

        Args:
            - model (libModel):
                The model loaded with libSBML.
            - gene_table (pd.DataFrame):
                The table with the genes to add. At least needs the columns
                *ncbiprotein* and *locus_tag*. Optional columns include
                *UniProt* amongst other.
        """

        # ncbiprotein | locus_tag | ...
        # work on a copy to ensure input stays the same
        gene_table = gene_table.copy()
        # gene_table.drop(columns=['ec-code'],inplace=True)

        # create gps from the table and add them to the model
        cols_for_refs = [
            _ for _ in REF_COL_GF_GENE_MAP.keys() if _ in gene_table.columns
        ]

        # create gp
        for idx, x in tqdm(
            gene_table.iterrows(), total=len(gene_table), desc="Adding genes to model"
        ):
            # get additional references
            references = dict()
            for dbname in cols_for_refs:
                references[dbname] = (x[dbname], True)
            create_gp(
                model, x["ncbiprotein"], locus_tag=x["locus_tag"], reference=references
            )
            self._statistics["genes"]["added"] += 1

    def add_gene_reac_associations_from_table(
        self, model: libModel, reac_table: pd.DataFrame
    ) -> None:
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
        reac_table = reac_table[["ncbiprotein", "add_to_GPR"]]
        reac_table = reac_table.explode("ncbiprotein")
        reac_table.drop_duplicates(subset="ncbiprotein", inplace=True)

        # add the genes to the corresponding GPRs
        for idx, row in reac_table.iterrows():
            # check, if G_+ncbiprotein in model
            # if yes, add gpr
            geneid = _f_gene_rev(row["ncbiprotein"])
            for reacid in row["add_to_GPR"]:
                current_reacid = "R_" + reacid
                current_mgids = [mgid for mgid in model_gene_ids if geneid in mgid]
                if current_mgids:
                    if len(current_mgids) == 1:
                        create_gpr(model.getReaction(current_reacid), current_mgids[0])
                    else:
                        mes = f"Found multiple matches for {geneid} in model: {current_mgids}. Belongs to reaction {current_reacid}."
                        warnings.warn(mes, UserWarning)
                # else, print warning
                else:
                    mes = f"Cannot find {geneid} in model. Should be added to {current_reacid}."
                    warnings.warn(mes, UserWarning)

    def add_reactions_from_table(
        self,
        model: cobra.Model,
        missing_reac_table: pd.DataFrame,
        formula_check: Literal["none", "existence", "wildcard", "strict"] = "existence",
        exclude_dna: bool = True,
        exclude_rna: bool = True,
        idprefix: str = "refineGEMs",
        namespace: Literal["BiGG"] = "BiGG",
    ) -> pd.DataFrame:
        """Helper function to add reactions to a model from the missing_reactions table
        (output of the chosen implementation of :py:meth:`~refinegems.classes.gapfill.GapFiller.find_missing_reactions`)

        Args:
            - model (cobra.Model):
                The model, loaded with COBRpy.
            - missing_reac_table (pd.DataFrame):
                The missing reactions table.
            - formula_check (Literal['none','existence','wildcard','strict'], optional):
                Param for checking metabolite formula before adding them to the model.
                For more information, refer to :py:func:`~refinegems.utility.entities.isreaction_complete`
                Defaults to 'existence'.
            - exclude_dna (bool, optional):
                Option to exclude reaction containing 'DNA' from being added to the model.
                Defaults to True.
            - exclude_rna (bool, optional):
                Option to exclude reaction containing 'RNA' from being added to the model.
                Defaults to True.
            - idprefix (str, optional):
                A prefix to use, if pseudo-IDs need to be created.
                Defaults to 'refineGEMs'.
            - namespace (Literal['BiGG'], optional):
                Namespace to use for the reactions and metabolites
                (and the model). Defaults to 'BiGG'.

        Raises:
            - TypeError: Unknown return type for reac param. Please contact the developers.

        Returns:
            pd.DataFrame:
                Table containing the information about which genes can now be added
                to reactions (use for GPR curation).
        """

        # reconstruct reactions
        # ---------------------
        for idx, row in tqdm(
            missing_reac_table.iterrows(),
            desc="Trying to add missing reacs",
            total=missing_reac_table.shape[0],
        ):
            # add EC number to references
            if row["reference"]:
                refs = row["reference"]
                if isinstance(refs, dict):
                    continue
                elif refs[0] == "{":
                    refs = refs.replace(r"'", r"\"")
                    refs = json.loads(refs)
                else:
                    refs = refs.split(r":")
                    refs = {refs[0]: refs[1]}
            else:
                refs = {}

            if row["ec-code"]:
                if "ec-code" in refs.keys():
                    if not isinstance(refs["ec-code"], list):
                        refs["ec-code"] = [refs["ec-code"]]
                    if isinstance(row["ec-code"], list):
                        refs["ec-code"] = list(set(refs["ec-code"] + row["ec-code"]))
                    elif isinstance(row["ec-code"], str):
                        refs["ec-code"] = list(
                            set(refs["ec-code"].append(row["ec-code"]))
                        )
                    else:
                        warnings.warn(
                            f'Unknown type for ec-code: {type(row["ec-code"])}',
                            UserWarning,
                        )
                else:
                    refs["ec-code"] = (
                        row["ec-code"]
                        if isinstance(row["ec-code"], list)
                        else [row["ec-code"]]
                    )

            # build reaction
            reac = None
            match row["via"]:
                # MetaNetX
                case "MetaNetX":
                    reac = build_reaction_mnx(
                        model,
                        row["id"],
                        reac_str=str(row["equation"]),
                        references=refs,
                        idprefix=idprefix,
                        namespace=namespace,
                    )
                # KEGG
                case "KEGG":
                    reac = build_reaction_kegg(
                        model,
                        row["id"],
                        reac_str=str(row["equation"]),
                        references=refs,
                        idprefix=idprefix,
                        namespace=namespace,
                    )
                # BiGG
                case "BiGG":
                    reac = build_reaction_bigg(
                        model,
                        row["id"],
                        references=refs,
                        idprefix=idprefix,
                        namespace=namespace,
                    )

                # Unknown database
                case _:
                    mes = f"""Unknown database name for reaction reconstruction: {row["via"]}\n
                    Reaction will not be reconstructed."""
                    warnings.warn(mes, UserWarning)

            # check output of reconstruction
            # ------------------------------
            # case 1: reconstruction was not possible
            if not reac:
                self._statistics["reactions"]["building failed"] += 1
                pass  # nothing to do here
            # case 2: reaction(s) found in model
            elif isinstance(reac, list):
                # add found names to the add_to_GPR column of the table
                current_gpr = missing_reac_table.loc[idx, "add_to_GPR"]
                if str(current_gpr) == "nan":
                    current_gpr = None
                if not current_gpr:
                    missing_reac_table.at[idx, "add_to_GPR"] = reac
                else:
                    missing_reac_table.at[idx, "add_to_GPR"] = list(
                        set(reac + current_gpr)
                    )
            # case 3: new reaction was generated
            elif isinstance(reac, cobra.Reaction):
                # validate reaction
                if isreaction_complete(
                    reac,
                    formula_check=formula_check,
                    exclude_dna=exclude_dna,
                    exclude_rna=exclude_rna,
                ):
                    # Extend reaction notes with information about the GapFiller
                    reac.notes["found with"] = f"refineGEMs GapFiller, {self._variety}"
                    # add reaction to model (if validation successful)
                    model.add_reactions([reac])
                    self._statistics["reactions"]["added"] += 1
                    # add reaction ID to table under add_to_GPR
                    current_gpr = missing_reac_table.loc[idx, "add_to_GPR"]
                    if pd.isna(current_gpr):
                        missing_reac_table.at[idx, "add_to_GPR"] = [reac.id]
                    else:
                        current_gpr.append(reac.id)
                        missing_reac_table.at[idx, "add_to_GPR"] = list(
                            set(current_gpr)
                        )

            # case 4: should never occur
            else:
                mes = f"Unknown return type for reac param. Please contact the developers."
                raise TypeError(mes)

        # save reactions, that could not be recontructed, for manual curation
        manual_curation_reacs = missing_reac_table[
            missing_reac_table["add_to_GPR"].isnull()
        ]
        self.manual_curation["reactions"]["building failed"] = manual_curation_reacs
        self._statistics["reactions"]["building failed"] = manual_curation_reacs[
            "id"
        ].nunique()
        # return the updated table with successfully reconstructed reaction ids
        # to enable adding the genes
        missing_gprs = missing_reac_table[~missing_reac_table["add_to_GPR"].isnull()]
        return missing_gprs

    def fill_model(self, model: Union[cobra.Model, libModel], **kwargs) -> libModel:
        """Based on a table of missing genes and missing reactions,
        fill the gaps in a model as good as possible automatically.

        .. note::

            This model rewrites and reloads the input model. Only the returned model
            has all the edits.

        Args:
            - model (Union[cobra.Model,libModel]):
                The model, either a libSBML or COBRApy model entity.
            - kwargs:
                Additional parameters to be passed to
                :py:meth:`~refinegems.classes.gapfill.GapFiller.add_reactions_from_table`.

        Raises:
            - TypeError: Unknown type of model.

        Returns:
            libModel:
                The gap-filled model.
        """
        # Step 0: Preparations
        # --------------------
        # load the correct type of model for the first step
        match model:
            case cobra.Model():
                with NamedTemporaryFile(suffix=".xml", delete=False) as tmp:
                    write_model_to_file(model, tmp.name)
                    model = load_model(tmp.name, "libsbml")
                os.remove(tmp.name)
            case libModel():
                pass
            case _:
                mes = f"Unknown type of model: {type(model)}"
                raise TypeError(mes)

        # Check if missing reactions found, if not return model
        if not isinstance(self.missing_reactions, pd.DataFrame) or self.missing_reactions.empty:
            return model

        # Filter out reactions without ncbiprotein
        self.manual_curation["reactions"]["no GPR"] = self.missing_reactions[
            self.missing_reactions["ncbiprotein"].isnull()
        ]
        self._statistics["reactions"]["unmappable"] += self.manual_curation[
            "reactions"
        ]["no GPR"]["id"].nunique()
        self.missing_reactions = self.missing_reactions[
            ~self.missing_reactions["ncbiprotein"].isnull()
        ]

        # Filter out genes without ncbiprotein
        self.manual_curation["genes"]["no ncbiprotein"] = self.missing_genes[
            self.missing_genes["ncbiprotein"].isnull()
        ]
        num_genes_no_ncbi_id = self.manual_curation["genes"]["no ncbiprotein"][
            "locus_tag"
        ].nunique()
        self._statistics["genes"]["unmappable"] += num_genes_no_ncbi_id
        self._statistics["genes"]["missing (mappable)"] -= num_genes_no_ncbi_id
        self.missing_genes = self.missing_genes[
            ~self.missing_genes["ncbiprotein"].isnull()
        ]

        # filter out duplicated genes to avoid duplicated IDs in the model
        if len(self.missing_genes) != len(self.missing_genes["ncbiprotein"].unique()):
            self.manual_curation["genes"]["duplicated (not added)"] = (
                self.missing_genes[
                    self.missing_genes.duplicated(subset=["ncbiprotein"])
                ]
            )
            self._statistics["genes"]["duplicated"] = self.manual_curation["genes"][
                "duplicated (not added)"
            ]["locus_tag"].nunique()
            self.missing_genes = self.missing_genes[
                ~self.missing_genes.duplicated(subset=["ncbiprotein"])
            ]

        # Step 1: Add genes to model whose reactions are already in it
        # -------------------------------------------------------------
        # filter the respective genes and reactions
        reacs_in_model = self.missing_reactions[
            ~(self.missing_reactions["add_to_GPR"].isnull())
        ]
        ncbiprot_with_reacs_in_model = [*chain(*list(reacs_in_model["ncbiprotein"]))]
        genes_with_reacs_in_model = self.missing_genes[
            self.missing_genes["ncbiprotein"].isin(ncbiprot_with_reacs_in_model)
        ]

        if len(genes_with_reacs_in_model) > 0:
            # add genes as gene products to model
            self.add_genes_from_table(model, genes_with_reacs_in_model)
            # extend gene production rules
            self.add_gene_reac_associations_from_table(model, reacs_in_model)

            # what remains:
            self.missing_reactions = self.missing_reactions[
                self.missing_reactions["add_to_GPR"].isnull()
            ]
            self.missing_genes = self.missing_genes[
                ~(self.missing_genes["ncbiprotein"].isin(ncbiprot_with_reacs_in_model))
            ]  #

        # Step 2: Add reactions to model, if reconstruction successful
        # ------------------------------------------------------------

        if len(self.missing_reactions) > 0:
            # re-load model with cobrapy
            with NamedTemporaryFile(suffix=".xml", delete=False) as tmp:
                write_model_to_file(model, tmp.name)
                model = load_model(tmp.name, "cobra")
            os.remove(tmp.name)

            # .......................
            # @DEBUG
            # if len(self.missing_reactions) > 10:
            #     self.missing_reactions = self.missing_reactions.sample(20)
            #     print('fill_model: Running in debugging mode')
            # .......................

            # add reactions to model
            missing_gprs = self.add_reactions_from_table(
                model, self.missing_reactions, **kwargs
            )

        # Step 3: Add GPRs + genes for the newly curated reactions
        # --------------------------------------------------------

        # re-load model with libsbml
        with NamedTemporaryFile(suffix=".xml", delete=False) as tmp:
            write_model_to_file(model, tmp.name)
            model = load_model(tmp.name, "libsbml")
        os.remove(tmp.name)

        try:
            if len(missing_gprs) > 0:
                # filter for genes for GPRs but not yet in model
                ncbiprot_with_reacs_in_model = [
                    *chain(*list(missing_gprs["ncbiprotein"]))
                ]
                genes_with_reacs_in_model = self.missing_genes[
                    self.missing_genes["ncbiprotein"].isin(ncbiprot_with_reacs_in_model)
                ]
                if len(genes_with_reacs_in_model) > 0:
                    # add genes as gene products to model
                    self.add_genes_from_table(model, genes_with_reacs_in_model)
                    # extend gene production rules
                    reacs_in_model = self.missing_reactions[
                        ~(self.missing_reactions["add_to_GPR"].isnull())
                    ]
                    self.add_gene_reac_associations_from_table(model, reacs_in_model)

                    self.missing_genes = self.missing_genes[
                        ~(
                            self.missing_genes["ncbiprotein"].isin(
                                ncbiprot_with_reacs_in_model
                            )
                        )
                    ]
        except NameError:
            warnings.warn("Something went wrong. Contact the developers (gapfillmodel)")

        # collect stats and stuff for manual curation
        self.manual_curation["genes"]["building failed"] = self.missing_genes
        self._statistics["genes"]["building failed"] = self.missing_genes[
            "locus_tag"
        ].nunique()

        # calculate remaining statistics
        self._statistics["reactions"]["missing (remaining)"] += self._statistics[
            "reactions"
        ]["building failed"]
        self._statistics["genes"]["missing (remaining)"] += self._statistics["genes"][
            "building failed"
        ]

        return model

    def report(
        self, dir: str, hide_zeros: bool = False, no_title: bool = False
    ) -> None:
        """Based on the previous gap-filling, save statistics and missing genes/reactions for manual curation.

        Args:
            - dir (str):
                Path to a directory to save the report to.
            - hide_zeros (bool, optional):
                Option to hide statistics with zero counts. Defaults to False.
            - with_title (bool, optional):
                Option to get figure without title. Defaults to False.
        """
        logging.warning(
            "Please keep in mind that all statistical values are determined while running the tool "
            + "and only unique values are tracked."
        )

        statistics_report = GapFillerReport(
            self._variety,
            deepcopy(self._statistics),
            self.manual_curation,
            hide_zeros,
            no_title,
        )
        statistics_report.save(Path(dir))


# --------------------
# Gapfilling with KEGG
# --------------------


class KEGGapFiller(GapFiller):
    """Based on a KEGG organism ID (corresponding to the organism of the model),
    find missing genes in the model and map them to reactions to try and fill the gaps
    found with the KEGG database.

    .. note::

        Please keep in mind that using this module requires a model containing the Genbank locus tags as labels.
        If your model does not conform to this you can use one of the functions
        :py:func:`~refinegems.curation.curate.polish_model` or
        :py:func:`~refinegems.curation.curate.extend_gp_annots_via_mapping_table`.

    .. hint::

        Due to the KEGG REST API this is relatively slow.

    Attributes:

    - GapFiller Attributes:
        All attributes of the parent class :py:class:`~refinegems.classes.gapfill.GapFiller`
    - organismid (str, required):
        Abbreviation of the organism in the KEGG database.

    """

    def __init__(self, organismid) -> None:
        super().__init__()
        self.organismid = organismid
        self._variety = "KEGG"

    def find_missing_genes(self, model: libModel):
        """Get the missing genes in model in comparison to the KEGG entry of the
        organism. Saves a table containing the missing genes
        to the attribute missing_genes.

        Format:

        ``orgid:locus | locus_tag | kegg.orthology | ec-code | ncbiprotein | uniprot``

        Args:
            - model (libModel):
                The model loaded with libSBML.

        """

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
                            if "kegg.genes" in uri:
                                genes_in_model.append(
                                    re.split(r"kegg.genes:|kegg.genes/", uri)[1]
                                )  # work with old/new pattern

            return pd.DataFrame(genes_in_model, columns=["orgid:locus"])

        # Step 1: get genes from model
        # ----------------------------
        genes_in_model = get_model_genes(model)

        # Step 2: get genes of organism from KEGG
        # ---------------------------------------
        gene_KEGG_list = KEGG().list(self.organismid)
        gene_KEGG_table = pd.read_table(io.StringIO(gene_KEGG_list), header=None)
        gene_KEGG_table.columns = ["orgid:locus", "CDS", "position", "protein"]
        self.full_gene_list = gene_KEGG_table
        gene_KEGG_table = gene_KEGG_table[["orgid:locus"]]

        # Statistics on full gene list based on KEGG
        self._statistics["genes"][f"total (based on {self._variety})"] = (
            self.full_gene_list["orgid:locus"].nunique()
            + int(self.full_gene_list["orgid:locus"].isna().sum())
        )

        # Step 3: KEGG vs. model genes -> get missing genes for model
        # ----------------------------
        genes_not_in_model = gene_KEGG_table[
            ~gene_KEGG_table["orgid:locus"].isin(genes_in_model["orgid:locus"])
        ]

        # Step 4: extract locus tag
        # -------------------------
        genes_not_in_model["locus_tag"] = (
            genes_not_in_model["orgid:locus"].str.split(r":").str[1]
        )

        # Step 5: map to EC via KEGG
        # --------------------------
        # @DEBUG .......................
        # genes_not_in_model = genes_not_in_model.iloc[0:50,:]
        # print(UserWarning('Running in debugging mode.'))
        # ..............................
        geneKEGG_mapping = pd.DataFrame.from_dict(
            list(genes_not_in_model["orgid:locus"].progress_apply(parse_KEGG_gene))
        )
        genes_not_in_model = genes_not_in_model.merge(
            geneKEGG_mapping, how="left", on="orgid:locus"
        )
        genes_not_in_model = genes_not_in_model.explode("ncbiprotein")

        # collect stats
        self._statistics["genes"]["missing (total)"] = genes_not_in_model[
            "locus_tag"
        ].nunique()

        self.missing_genes = genes_not_in_model

    def find_missing_reactions(self, model: cobra.Model, threshold_add_reacs: int = 5):

        # Step 1: filter missing gene list + extract ECs
        # ----------------------------------------------
        # Filter out unmappable genes that have whether NCBI protein IDs nor EC numbers
        mask = (
            self.missing_genes["ec-code"].isnull()
            & self.missing_genes["ncbiprotein"].isnull()
        )
        self.manual_curation["genes"]["no ncbiprotein, no EC"] = self.missing_genes[
            mask
        ]
        self.missing_genes = self.missing_genes[~mask]
        self._statistics["genes"]["unmappable"] = self.manual_curation["genes"][
            "no ncbiprotein, no EC"
        ]["locus_tag"].nunique()

        # Filter out remaining unmappable genes due to missing EC numbers
        mask = self.missing_genes["ec-code"].isnull()
        self.manual_curation["genes"]["no EC"] = self.missing_genes[mask]
        self._statistics["genes"]["unmappable"] += self.manual_curation["genes"][
            "no EC"
        ]["locus_tag"].nunique()
        self.missing_genes = self.missing_genes[~mask]

        # collect overall remaining statistics on genes
        self._statistics["genes"]["missing (mappable)"] = self.missing_genes[
            "locus_tag"
        ].nunique()
        self._statistics["genes"]["missing (remaining)"] = self._statistics["genes"][
            "unmappable"
        ]

        # get relevant infos for reacs
        self.missing_reactions = self.missing_genes[["ec-code", "ncbiprotein"]]

        # check, if any automatic gapfilling is possible
        if self.missing_reactions.empty:
            logging.warning(
                f"No missing reactions for the provided model {model.id} were found via {self._variety}."
            )
            return None
        # transform table into EC-number vs. list of NCBI protein IDs
        eccode = (
            self.missing_reactions["ec-code"]
            .apply(pd.Series)
            .reset_index()
            .melt(id_vars="index")
            .dropna()[["index", "value"]]
            .set_index("index")
        )
        ncbiprot = (
            self.missing_reactions["ncbiprotein"]
            .apply(pd.Series)
            .reset_index()
            .melt(id_vars="index")
            .dropna()[["index", "value"]]
            .set_index("index")
        )
        self.missing_reactions = pd.merge(
            eccode, ncbiprot, left_index=True, right_index=True
        ).rename(columns={"value_x": "ec-code", "value_y": "ncbiprotein"})
        self.missing_reactions = (
            self.missing_reactions.groupby(self.missing_reactions["ec-code"])
            .aggregate({"ncbiprotein": "unique"})
            .reset_index()
        )

        # @DEBUG .......................
        # self.missing_reactions = self.missing_reactions.iloc[10:30,:]
        # print(UserWarning('Running in debugging mode.'))
        # ..............................

        # Step 2: map EC to reaction(s) if possible
        # -----------------------------------------
        # via MNX, BiGG, KEGG
        reacs_mapped = map_ec_to_reac(self.missing_reactions, threshold_add_reacs)

        # Statistics on missing based on genes
        self._statistics["reactions"]["missing (based on genes)"] = reacs_mapped[
            "id"
        ].nunique() + int(reacs_mapped["id"].isna().sum())

        # Step 3: clean and map to model reactions
        # ----------------------------------------
        # need manual curation
        self.manual_curation["reactions"]["no ID"] = reacs_mapped[
            reacs_mapped["id"].isnull()
        ]
        self._statistics["reactions"]["unmappable"] = int(
            self.manual_curation["reactions"]["no ID"]["id"].isna().sum()
        )

        # map to model reactions
        reacs_mapped = reacs_mapped[~reacs_mapped["id"].isnull()]
        gpr = reacs_mapped.apply(
            lambda x: self._find_reac_in_model(model, x["ec-code"], x["id"], x["via"]),
            axis=1,
        )
        if type(gpr) == object:
            reacs_mapped["add_to_GPR"] = gpr
        else:
            reacs_mapped["add_to_GPR"] = None

        # statistics on reactions already in model
        self._statistics["reactions"]["already in model"] = reacs_mapped[
            ~reacs_mapped["add_to_GPR"].isnull()
        ]["id"].nunique()

        # statistics on mapped reactions
        self._statistics["reactions"]["mapped to MNX"] = reacs_mapped[
            reacs_mapped["add_to_GPR"].isnull()
        ][reacs_mapped["via"] == "MetaNetX"]["id"].nunique()
        self._statistics["reactions"]["mapped to BiGG"] = reacs_mapped[
            reacs_mapped["add_to_GPR"].isnull()
        ][reacs_mapped["via"] == "BiGG"]["id"].nunique()
        self._statistics["reactions"]["mapped to KEGG"] = reacs_mapped[
            reacs_mapped["add_to_GPR"].isnull()
        ][reacs_mapped["via"] == "KEGG"]["id"].nunique()

        # calculate remaining statistics
        self._statistics["reactions"]["missing (total)"] = (
            self._statistics["reactions"]["missing (based on genes)"]
            - self._statistics["reactions"]["already in model"]
        )
        self._statistics["reactions"]["missing (remaining)"] = self._statistics[
            "reactions"
        ]["unmappable"]

        # return missing_reactions
        self.missing_reactions = reacs_mapped


# ----------------------
# Gapfilling with BioCyc
# ----------------------
class BioCycGapFiller(GapFiller):
    """
    | Based on a SmartTable with information on the genes and a SmartTable with
    | information on the reactions of the organism of the model, this class
    | finds missing genes in the model and maps them to reactions to try and
    | fill the gaps found with the BioCyc gene SmartTable.
    |
    | For specifications on the SmartTables see the attributes `biocyc_gene_tbl`
    | & `biocyc_reacs_tbl`

    .. note::

        Please keep in mind that using this module requires a model containing the Genbank locus tags as labels.
        If your model does not conform to this you can use one of the functions
        :py:func:`~refinegems.curation.curate.polish_model` or
        :py:func:`~refinegems.curation.curate.extend_gp_annots_via_mapping_table`.

    Attributes:
        - GapFiller Attributes:
            All attributes of the parent class :py:class:`~refinegems.classes.gapfill.GapFiller`
        - biocyc_gene_tbl_path (str, required):
            Path to organism-specific SmartTable for genes from BioCyc;
            Should contain the columns: ``Accession-2 | Reactions of gene``
        - biocyc_reacs_tbl_path (str, required):
            Path to organism-specific SmartTable for reactions from BioCyc;
            Should contain the columns:
            ``Reaction | Object ID | EC-Number | Spontaneous?``
        - gff (str, required):
            Path to organism-specific GFF file
    """

    def __init__(
        self, biocyc_gene_tbl_path: str, biocyc_reacs_tbl_path: str, gff: str
    ) -> None:
        super().__init__()
        self.full_gene_list = biocyc_gene_tbl_path
        self.biocyc_rxn_tbl = biocyc_reacs_tbl_path
        self._gff = gff
        self._variety = "BioCyc"

    @property
    def full_gene_list(self):
        """
        | Get or set the current BioCyc Gene table.
        | While setting the provided path for a TSV file from BioCyc with the
        | columns ``'Accession-2' | 'Reactions of gene'`` is parsed and the
        | content of the file is stored in a DataFrame containing all rows where
        | a 'Reactions of gene' exists.

        .. hint::

            Please keep in mind that the column Accession-2 needs to contain Genbank locus tags. If that is not the case
            for your organism use the correct column from BioCyc and rename it accordingly.

        """
        return self._full_gene_list

    @full_gene_list.setter
    def full_gene_list(self, biocyc_gene_tbl_path: str):
        # Read table
        biocyc_genes = pd.read_table(
            biocyc_gene_tbl_path,
            usecols=["Accession-2", "Reactions of gene"],
            dtype=str,
        )

        # Rename columns for further use
        biocyc_genes.rename(
            columns={"Accession-2": "locus_tag", "Reactions of gene": "id"},
            inplace=True,
        )

        # Turn empty strings into None
        biocyc_genes.replace(r"", None, inplace=True)

        # Drop only complete empty rows
        biocyc_genes.dropna(how="all", inplace=True)

        # Statistics on full gene list based on BioCyc
        self._statistics["genes"][f"total (based on {self._variety})"] = biocyc_genes[
            "locus_tag"
        ].nunique() + int(self.missing_genes["locus_tag"].isna().sum())

        self.full_gene_list = biocyc_genes

    @property
    def biocyc_rxn_tbl(self):
        """
        | Get or set the current BioCyc Reaction table.
        | While setting the provided path for a TSV file from BioCyc with the
        | columns ``'Reaction' | 'Object ID' | 'EC-Number' | 'Spontaneous?'`` is
        | parsed and the content of the file is stored in a DataFrame.
        """
        return self._biocyc_rxn_tbl

    @biocyc_rxn_tbl.setter
    def biocyc_rxn_tbl(self, biocyc_reacs_tbl_path: str) -> pd.DataFrame:
        # Read table
        biocyc_rxns = pd.read_table(
            biocyc_reacs_tbl_path,
            usecols=["Reaction", "Object ID", "EC-Number", "Spontaneous?"],
            dtype=str,
        )

        # Rename columns for further use
        biocyc_rxns.rename(
            columns={
                "Reaction": "equation",
                "Object ID": "id",
                "EC-Number": "ec-code",
                "Spontaneous?": "is_spontaneous",
            },
            inplace=True,
        )

        # Turn empty strings into NaNs
        biocyc_rxns.replace(r"", None, inplace=True)

        # Set entries in is_spontaneous to booleans &
        # specify empty entries in 'is_spontaneous' as False
        biocyc_rxns["is_spontaneous"].replace({r"T": True, r"F": False}, inplace=True)
        biocyc_rxns["is_spontaneous"] = biocyc_rxns["is_spontaneous"].fillna(False)

        self._biocyc_rxn_tbl = biocyc_rxns

    def find_missing_genes(self, model: libModel):
        """Retrieves the missing genes and reactions from the BioCyc table
        according to the 'Accession-2' identifiers (locus_tags)

        Args:
            - model (libModel):
                Model loaded with libSBML
        """

        # Step 1: get genes from model
        # ----------------------------
        geneps_in_model = [
            _.getLabel() for _ in model.getPlugin(0).getListOfGeneProducts()
        ]

        # Step 2: Get genes of organism from BioCyc
        # -----------------------------------------
        # See self._biocyc_gene_tbl

        # Step 3: BioCyc vs. model genes -> get missing genes for model
        # -------------------------------------------------------------
        self.missing_genes = self.biocyc_gene_tbl[
            ~self.biocyc_gene_tbl["locus_tag"].isin(geneps_in_model)
        ]

        # Step 4: Get amount of missing genes in total
        # --------------------------------------------
        self._statistics["genes"]["missing (total)"] = self.missing_genes[
            "locus_tag"
        ].nunique()

        # Step 5: Filter results
        # ----------------------
        # Save not mappable genes due to no reaction ID
        self.manual_curation["genes"]["no reaction ID"] = self.missing_genes[
            self.missing_genes["id"].isnull()
        ]

        # Add amount of unmappable genes to statistics
        self._statistics["genes"]["unmappable"] = self.manual_curation["genes"][
            "no reaction ID"
        ]["locus_tag"].nunique()

        # Remove all rows where 'id' is None
        self.missing_genes.dropna(subset="id", inplace=True)

        # Step 6: Get ncbiprotein IDs
        # ---------------------------
        # Parse GFF file to obtain locus_tag2ncbiportein mapping for all CDS
        locus_tag2ncbiprotein_df = parse_gff_for_cds(
            self._gff,
            {"locus_tag": "locus_tag", "protein_id": "ncbiprotein", "product": "name"},
        )
        locus_tag2ncbiprotein_df = locus_tag2ncbiprotein_df.explode("ncbiprotein")
        locus_tag2ncbiprotein_df = locus_tag2ncbiprotein_df.explode("name")

        # Get the complete missing genes dataframe with the ncbiprotein IDs
        self.missing_genes = self.missing_genes.merge(
            locus_tag2ncbiprotein_df, on="locus_tag"
        )

        # Step 7: Get amount of missing genes in total
        # --------------------------------------------
        self._statistics["genes"]["missing (mappable)"] = self.missing_genes[
            "locus_tag"
        ].nunique()
        self._statistics["genes"]["missing (remaining)"] = self._statistics["genes"][
            "unmappable"
        ]

    def find_missing_reactions(self, model: cobra.Model):
        """Retrieves the missing reactions with more information like the
        equation, EC code, etc. according to the missing genes

        Args:
            - model (cobra.Model):
                Model loaded with COBRApy
        """

        # Step 1: filter missing gene list + extract ECs
        # ----------------------------------------------
        # Drop locus tag column as not needed for here
        missing_genes = self.missing_genes.drop(["locus_tag", "name"], axis=1)

        # Expand missing genes result table to merge with Biocyc reactions table
        missing_genes = pd.DataFrame(
            missing_genes["id"].str.split(r"//").tolist(),
            index=missing_genes["ncbiprotein"],
        ).stack()
        missing_genes = missing_genes.reset_index([0, "ncbiprotein"])
        missing_genes.columns = ["ncbiprotein", "id"]
        missing_genes["id"] = missing_genes["id"].str.strip()

        # Turn ncbiprotein column into lists of ncbiprotein IDs per reaction
        ncbiprotein_as_list = (
            missing_genes.groupby("id")["ncbiprotein"]
            .apply(list)
            .reset_index(name="ncbiprotein")
        )
        missing_genes.drop("ncbiprotein", axis=1, inplace=True)
        missing_genes = ncbiprotein_as_list.merge(missing_genes, on="id")
        missing_genes["ncbiprotein"] = missing_genes["ncbiprotein"].apply(
            lambda x: x if not None in x else list(filter(None, x))
        )

        # Drop duplicates to get the unique BioCyc IDs
        missing_genes.drop_duplicates(subset="id", inplace=True)

        # Get missing reactions from missing genes
        self.missing_reactions = missing_genes.merge(self.biocyc_rxn_tbl, on="id")

        # Turn ec-code entries with '//' into lists, remove prefix 'EC-' & get unique ec-codes
        self.missing_reactions["ec-code"] = (
            self.missing_reactions["ec-code"]
            .str.replace(r"EC-", r"")
            .str.split(r"\s*//\s*")
        )
        self.missing_reactions["ec-code"] = (
            self.missing_reactions["ec-code"]
            .fillna({i: [] for i in self.missing_reactions.index})
            .map(set)
            .map(list)
        )

        # Add 'G_spontaneous' as gene product if marked as spontaneous &
        # drop is_spontaneous column
        self.missing_reactions["ncbiprotein"] = self.missing_reactions.apply(
            lambda x: (
                x["ncbiprotein"]
                if not x["is_spontaneous"]
                else x["ncbiprotein"].append("spontaneous")
            ),
            axis=1,
        )
        self.missing_reactions.drop("is_spontaneous", axis=1, inplace=True)

        # Step 2: Get amount of missing reactions from BioCyc for statistics
        # ------------------------------------------------------------------
        self._statistics["reactions"]["missing (based on genes)"] = (
            self.missing_reactions["id"].nunique()
        )

        # Step 3: Map BioCyc to model reactions & cleanup
        # -----------------------------------------------
        # Add column 'via'
        self.missing_reactions["via"] = "BioCyc"

        # Filter reacs for already in model
        self.missing_reactions["add_to_GPR"] = self.missing_reactions.apply(
            lambda x: self._find_reac_in_model(model, x["ec-code"], x["id"], x["via"]),
            axis=1,
        )

        # Add column 'reference'
        self.missing_reactions["reference"] = None

        # Drop rows if 'id' is NaN
        self.manual_curation["reactions"]["no ID"] = self.missing_reactions[
            self.missing_reactions["id"].isnull()
        ]

        self._statistics["reactions"]["unmappable"] = int(
            self.manual_curation["reactions"]["no ID"]["id"].isna().sum()
        )
        self.missing_reactions = self.missing_reactions[
            ~self.missing_reactions["id"].isnull()
        ]

        # check, if any automatic gapfilling is possible
        if self.missing_reactions.empty:
            logging.warning(
                f"No missing reactions for the provided model {model.id} were found via {self._variety}."
            )
            return None

        # Step 4: Map missing reactions without entries in column 'add_to_GPR'
        #         to other databases to get a parsable reaction equation
        # --------------------------------------------------------------------
        # Map to MetaNetX/BiGG
        mapped_reacs = map_biocyc_to_reac(self.missing_reactions)

        # Filter reacs for already in model
        mapped_reacs["add_to_GPR"] = mapped_reacs.apply(
            lambda x: self._find_reac_in_model(model, x["ec-code"], x["id"], x["via"]),
            axis=1,
        )

        # Step 5: Get results
        # -------------------
        # statistics on reactions already in model
        self._statistics["reactions"]["already in model"] = mapped_reacs[
            ~mapped_reacs["add_to_GPR"].isnull()
        ]["id"].nunique()

        # statistics on mapped reactions
        self._statistics["reactions"]["mapped to MNX"] = mapped_reacs[
            mapped_reacs["via"] == "MetaNetX"
        ]["id"].nunique()
        self._statistics["reactions"]["mapped to BiGG"] = mapped_reacs[
            mapped_reacs["via"] == "BiGG"
        ]["id"].nunique()
        self._statistics["reactions"]["mapped to KEGG"] = mapped_reacs[
            mapped_reacs["via"] == "KEGG"
        ]["id"].nunique()

        # Split missing reactios based on entries in 'via' & 'add_to_GPR'
        mask = (mapped_reacs["via"] == "BioCyc") & (mapped_reacs["add_to_GPR"].isnull())

        # DataFrame with unmappable BioCyc IDs & No entries in 'add_to_GPR'
        self.manual_curation["reactions"]["no mapping"] = mapped_reacs[mask]
        self._statistics["reactions"]["unmappable"] += self.manual_curation[
            "reactions"
        ]["no mapping"]["id"].nunique()

        # calculate remaining statistics
        self._statistics["reactions"]["missing (total)"] = (
            self._statistics["reactions"]["missing (based on genes)"]
            - self._statistics["reactions"]["already in model"]
        )
        self._statistics["reactions"]["missing (remaining)"] = self._statistics[
            "reactions"
        ]["unmappable"]

        # Mapped reactions
        self.missing_reactions = mapped_reacs[~mask]


# -------------------------------------
# GapFilling with GFF and DMND database
# -------------------------------------


class GeneGapFiller(GapFiller):
    """Find gaps in the model using the GFF file of the underlying genome
    and a DMND database and optionally NCBI.

    This gap filling approach tries to identify missing genes from the GFF file
    and uses DIAMOND to run a blastp search for homologs against the DMND database

    .. note::

        Please keep in mind that using this module requires a model containing the Genbank locus tags as labels.
        If your model does not conform to this you can use one of the functions
        :py:func:`~refinegems.curation.curate.polish_model` or
        :py:func:`~refinegems.curation.curate.extend_gp_annots_via_mapping_table`.

    .. hint::

        Files required for the swissprot approach can be downloaded with
        :py:func:`~refinegems.utility.set_up.download_url`

    Attributes:
        - GapFiller Attributes:
            All attributes of the parent class :py:class:`~refinegems.classes.gapfill.GapFiller`
    """

    GFF_COLS = {
        "locus_tag": "locus_tag",
        "eC_number": "ec-code",
        "protein_id": "ncbiprotein",
    }  # :meta:

    def __init__(self) -> None:
        super().__init__()
        self._variety = "GFF"

    def find_missing_genes(self, gffpath: Union[str, Path], model: libModel):
        """Find missing genes by comparing the CDS regions written in the GFF
        with the GeneProduct entities in the model.

        Args:
            - gffpath (Union[str, Path]):
                Path to a GFF file (corresponding to the model).
            - model (libModel):
                The model loaded with libSBML.
        """

        # get all CDS from gff
        self.full_gene_list = parse_gff_for_cds(gffpath, self.GFF_COLS)
        # Statistics on full gene list based on GFF
        self._statistics["genes"][f"total (based on {self._variety})"] = (
            self.full_gene_list["locus_tag"].nunique()
            + int(self.full_gene_list["locus_tag"].isna().sum())
        )
        # get all genes from model by locus tag
        model_locustags = [
            _f_gene(g.getLabel()) for g in model.getPlugin(0).getListOfGeneProducts()
        ]
        # filter
        self.missing_genes = self.full_gene_list.loc[
            ~self.full_gene_list["locus_tag"].isin(model_locustags)
        ]
        # formatting
        for col in self.GFF_COLS.values():
            if col not in self.missing_genes.columns:
                self.missing_genes[col] = None

        # collect stats
        self._statistics["genes"]["missing (total)"] = self.missing_genes[
            "locus_tag"
        ].nunique() + int(self.missing_genes["locus_tag"].isna().sum())

        # save genes with no locus tag for manual curation
        if "ncbiprotein" in self.missing_genes.columns:
            self.manual_curation["genes"]["gff no locus tag"] = self.missing_genes[
                self.missing_genes["locus_tag"].isnull()
            ]["ncbiprotein"]
        else:
            self.missing_genes["ncbiprotein"] = None

        # Statistics on unmappable genes
        self._statistics["genes"]["unmappable"] = int(
            self.missing_genes["locus_tag"].isna().sum()
        )  # no locus tag

        # formatting
        # ncbiprotein | locus_tag | ec-code
        self.missing_genes = self.missing_genes[
            ~self.missing_genes["locus_tag"].isnull()
        ]
        self.missing_genes = self.missing_genes.explode("ncbiprotein")

    def find_missing_reactions(
        self,
        model: cobra.Model,
        # prefix for pseudo ncbiprotein ids
        prefix: str = "refinegems",
        type_db: Literal["swissprot", "user"] = "swissprot",
        # database information
        fasta: str = None,
        dmnd_db: str = None,
        map_db: str = None,
        # SwissProt - NCBI params
        mail: str = None,
        check_NCBI: bool = False,
        # other params
        threshold_add_reacs: int = 5,
        # further optional params for the mapping
        **kwargs,
    ) -> None:
        """Find missing reactions in the model by blasting the missing genes
        against the SwissProt database and mapping the results to EC/BRENDA.

        Optionally, query the protein accession numbers against NCBI
        (increases runtime significantly).

        .. hint::

            For more information on how to get the SwissProt files, please see
            :py:func:`~refinegems.utility.set_up.download_url`.

        Args:
            - model (cobra.Model):
                The model loaded with COBRApy.
            - prefix (str, optional):
                Prefix for gene IDs in the model, the have to be generated
                randomly, as no ID from the chosen namespace (usually NCBI protein)
                has been found. Defaults to 'refinegems'.
            - mail (str, optional):
                Mail address for the query against NCBI.
                If not set, skips the mapping.
                Defaults to None.
            - check_NCBI (bool, optional):
                If set to True, checking the gene IDs / NCBI protein IDs against the
                NCBI database is enabled. Else, this step is skipped to reduce runtime.
                Only usable with SwissProt as database.
                Defaults to False.
            - fasta (str, optional):
                The protein FASTA file of the organism the model was build on.
                Required for the searchh against SwissProt.
                Defaults to None.
            - type_db (Literal['swissprot','user'], optional):
                Database to search against.
                Choose 'swissprot' for SwissProt or 'user' for a user defined database.
                Defaults to 'swissprot'.
            - dmnd_db (str, optional):
                Path to the DIAMOND database containing the protein sequences of SwissProt.
                Required for the search against SwissProt or rhe users own DIAMOND dn.
                Defaults to None.
            - map_db (str, optional):
                Path to the SwissProt mapping file.
                Required for the search against SwissProt.
                Greatly decreases runtime for running the DIAMOND search.

                ..note::
                    The mapping depends on the chosen database.

                Defaults to None.
            - threshold_add_reacs (int, optional):
                Threshold for the amount of reactions to add to the model.
                Defaults to 5.
            - **kwargs:
                Further optional parameters for the mapping,
                e.g. outdir, sens, cov, t, pid, etc.
                For more information see :py:func:`refinegems.utility.db_access.map_to_homologs`
                in case of type_db = 'user' or :py:func:`refinegems.utility.db_access.get_ec_via_swissprot`
                in case of type_db = 'swissprot'.
        """
        self._variety = f"{self._variety} + {type_db}"

        # try to identfy missing ECs
        # --------------------------
        case_1 = self.missing_genes[self.missing_genes["ec-code"].isnull()]
        not_case_1 = self.missing_genes[~self.missing_genes["ec-code"].isnull()]
        if len(case_1) > 0:

            # BLAST against a dmnd to retrieve EC numbers
            # +++++++++++++++++++++++++++++++++++++++++++

            match type_db:
                # type_db = swissprot: BLAST against SwissProt
                # -> BLAST (DIAMOND) against SwissProt to get EC/BRENDA
                case "swissprot":
                    if fasta and dmnd_db and map_db:
                        case_1_mapped = get_ec_via_swissprot(
                            fasta, dmnd_db, case_1, map_db, **kwargs
                        )  # further optional params for the mapping
                        case_1.drop("ec-code", inplace=True, axis=1)
                        case_1 = case_1.merge(case_1_mapped, on="locus_tag", how="left")
                        not_case_1["UniProt"] = None
                        # still no EC but ncbiprotein
                        #       -> access ncbi for ec (optional)
                        # @DEBUG .......................
                        # mapped_reacs = mapped_reacs.iloc[300:350,:]
                        # print(UserWarning('Running in debugging mode.'))
                        # ..............................
                        if check_NCBI and mail:
                            mapped_reacs["ec-code"] = mapped_reacs.progress_apply(
                                lambda x: (
                                    get_ec_from_ncbi(mail, x["ncbiprotein"])
                                    if not x["ec-code"]
                                    and not x["ncbiprotein"].isnull()
                                    else x["ec-code"]
                                ),
                                axis=1,
                            )

                # type_db = user: BLAST against user defined database
                case "user":
                    if fasta and dmnd_db:
                        case_1 = map_to_homologs(
                            fasta, dmnd_db, case_1, map_db, email=mail, **kwargs
                        )

                case _:
                    raise ValueError(
                        "Please choose one of the valid database types 'swissprot' or 'user' for the GeneGapFiller."
                    )

        # no ncbiprotein, no EC
        self.manual_curation["genes"]["no ncbiprotein, no EC"] = case_1[
            case_1["ncbiprotein"].isnull() & case_1["ec-code"].isnull()
        ]
        self._statistics["genes"]["unmappable"] += self.manual_curation["genes"][
            "no ncbiprotein, no EC"
        ]["locus_tag"].nunique()

        # combine with the rest
        mapped_reacs = pd.concat(
            [
                case_1[~(case_1["ncbiprotein"].isnull() & case_1["ec-code"].isnull())],
                not_case_1,
            ]
        )

        # convert NaNs to None
        mapped_reacs.mask(mapped_reacs.isnull(), other=None, inplace=True)

        # save entries with no EC for manual curation
        self.manual_curation["genes"]["no EC"] = mapped_reacs[
            mapped_reacs["ec-code"].isnull()
        ]
        self._statistics["genes"]["unmappable"] += self.manual_curation["genes"][
            "no EC"
        ]["locus_tag"].nunique()
        mapped_reacs = mapped_reacs[~mapped_reacs["ec-code"].isnull()]

        # check, if any automatic gapfilling is still possible
        if mapped_reacs.empty:
            logging.warning(
                f"No missing reactions for the provided model {model.id} were found via {self._variety}."
            )
            return None

        # create pseudoids for entries with no ncbiprotein id
        mapped_reacs["ncbiprotein"] = mapped_reacs.apply(
            lambda x: f'{x["locus_tag"]}' if not x["ncbiprotein"] else x["ncbiprotein"],
            axis=1,
        )

        # EC found
        # ---------
        # update the gene information
        updated_missing_genes = mapped_reacs.copy()

        # reformat missing reacs
        if type_db == "swissprot":
            mapped_reacs.drop(["locus_tag"], inplace=True, axis=1)

        # transform table into EC-number vs. list of NCBI protein IDs
        eccode = (
            mapped_reacs["ec-code"]
            .apply(pd.Series)
            .reset_index()
            .melt(id_vars="index")
            .dropna()[["index", "value"]]
            .set_index("index")
        )
        ncbiprot = (
            mapped_reacs["ncbiprotein"]
            .apply(pd.Series)
            .reset_index()
            .melt(id_vars="index")
            .dropna()[["index", "value"]]
            .set_index("index")
        )
        mapped_reacs = pd.merge(
            eccode, ncbiprot, left_index=True, right_index=True
        ).rename(columns={"value_x": "ec-code", "value_y": "ncbiprotein"})
        mapped_reacs = (
            mapped_reacs.groupby(mapped_reacs["ec-code"])
            .aggregate({"ncbiprotein": "unique"})
            .reset_index()
        )

        # map EC to reactions
        mapped_reacs = map_ec_to_reac(
            mapped_reacs[["ec-code", "ncbiprotein"]], threshold_add_reacs
        )

        self._statistics["reactions"]["missing (based on genes)"] = mapped_reacs[
            "id"
        ].nunique()

        # save for manual curation
        self.manual_curation["reactions"]["no mapping"] = mapped_reacs[
            mapped_reacs["id"].isnull()
        ]
        self._statistics["reactions"]["unmappable"] = int(
            self.manual_curation["reactions"]["no mapping"]["id"].isna().sum()
        )
        # map to model
        mapped_reacs = mapped_reacs[~mapped_reacs["id"].isnull()]
        mapped_reacs["add_to_GPR"] = mapped_reacs.apply(
            lambda x: self._find_reac_in_model(model, x["ec-code"], x["id"], x["via"]),
            axis=1,
        )

        # Track all remaining statistics
        # statistics on reactions already in model
        self._statistics["reactions"]["already in model"] = mapped_reacs[
            ~mapped_reacs["add_to_GPR"].isnull()
        ]["id"].nunique()

        # statistics on mapped reactions
        self._statistics["reactions"]["mapped to MNX"] = mapped_reacs[
            mapped_reacs["add_to_GPR"].isnull()
        ][mapped_reacs["via"] == "MetaNetX"]["id"].nunique()
        self._statistics["reactions"]["mapped to BiGG"] = mapped_reacs[
            mapped_reacs["add_to_GPR"].isnull()
        ][mapped_reacs["via"] == "BiGG"]["id"].nunique()
        self._statistics["reactions"]["mapped to KEGG"] = mapped_reacs[
            mapped_reacs["add_to_GPR"].isnull()
        ][mapped_reacs["via"] == "KEGG"]["id"].nunique()

        # calculate remaining statistics
        self._statistics["reactions"]["missing (total)"] = (
            self._statistics["reactions"]["missing (based on genes)"]
            - self._statistics["reactions"]["already in model"]
        )
        self._statistics["reactions"]["missing (remaining)"] = self._statistics[
            "reactions"
        ]["unmappable"]
        self._statistics["genes"]["missing (mappable)"] = updated_missing_genes[
            "locus_tag"
        ].nunique()
        self._statistics["genes"]["missing (remaining)"] = self._statistics["genes"][
            "unmappable"
        ]

        # update attributes
        if type_db == "user":
            updated_missing_genes = updated_missing_genes.drop(
                columns=["locus_tag_ref", "old_locus_tag"], axis=1
            )
        self.missing_genes = updated_missing_genes
        self.missing_reactions = mapped_reacs


############################################################################
# functions (2)
############################################################################


# Gap-filling via medium/media
# ----------------------------
def single_cobra_gapfill(
    model: cobra.Model,
    universal: cobra.Model,
    medium: Medium,
    namespace: Literal["BiGG"] = "BiGG",
    growth_threshold: float = 0.05,
) -> Union[list[str], bool]:
    """Attempt gapfilling (with COBRApy) for a given model to allow growth on a given
    medium.

    Args:
        - model (cobra.Model):
            The model to perform gapfilling on.
        - universal (cobra.Model):
            A model with reactions to be potentially used for the gapfilling.
        - medium (Medium):
            A medium the model should grow on.
        - namespace (Literal['BiGG'], optional): Namespace to use for the model.
            Defaults to 'BiGG'.
        - growth_threshold (float, optional):  Minimal rate for the model to be considered growing.
            Defaults to 0.05.

    Returns:
        Union[list[str],True]:
            List of reactions to be added to the model to allow growth
            or True, if the model already grows.
    """

    # perform the gapfilling
    solution = []
    with model as model_copy:
        # set medium model should grow on
        medium_to_model(model_copy, medium, namespace, double_o2=False)
        # if model does not show growth (depending on threshold), perform gapfilling
        if model_copy.optimize().objective_value < growth_threshold:
            try:
                solution = cobra.flux_analysis.gapfill(
                    model_copy,
                    universal,
                    lower_bound=growth_threshold,
                    demand_reactions=False,
                )
            except cobra.exceptions.Infeasible:
                warnings.warn(
                    f"Gapfilling for medium {medium.name} failed. Manual curation required."
                )
        else:
            print(
                f"Model already grows on medium {medium.name} with objective value of {model_copy.optimize().objective_value}"
            )
            return True

    return solution


def cobra_gapfill_wrapper(
    model: cobra.Model,
    universal: cobra.Model,
    medium: Medium,
    namespace: Literal["BiGG"] = "BiGG",
    growth_threshold: float = 0.05,
    iterations: int = 3,
    chunk_size: int = 10000,
) -> cobra.Model:
    """Wrapper for :py:func:`~refinegems.classes.gapfill.single_cobra_gapfill`.

    Either use the full set of reactions in universal model by setting iteration to
    0 or None or use them in randomized chunks for faster runtime (useful
    on laptops). Note: when using the second option, be aware that this does not
    test all reaction combinations exhaustively (heuristic approach!!!).

    Args:
        - model (cobra.Model):
            The model to perform gapfilling on.
        - universal (cobra.Model):
            A model with reactions to be potentially used for the gapfilling.
        - medium (Medium):
            A medium the model should grow on.
        - namespace (Literal['BiGG'], optional):
            Namespace to use for the model.
            Options include 'BiGG'.
            Defaults to 'BiGG'.
        - growth_threshold (float, optional):
            Growth threshold for the gapfilling.
            Defaults to 0.05.
        - iterations (int, optional):
            Number of iterations for the heuristic version of the gapfilling.
            If 0 or None is given, uses full set of reactions.
            Defaults to 3.
        - chunk_size (int, optional):
            Number of reactions to be used for gapfilling at the same time.
            If None or 0 is given, use full set, not heuristic.
            Defaults to 10000.

    Returns:
        cobra.Model:
            The gapfilled model, if a solution was found.
    """

    solution = []

    # run a heuristic approach:
    #     for a given number of iterations, use a subset (chunk_size) of
    #     reactions for the gapfilling
    if (iterations and iterations > 0) and (chunk_size and chunk_size > 0):

        num_reac = len(model.reactions)
        # for each iteration
        for i in range(iterations):
            not_used = [_ for _ in range(0, num_reac)]

            # divide reactions in random subsets
            for chunk in range(math.ceil(num_reac / chunk_size)):
                if len(not_used) > chunk_size:
                    rng = np.random.default_rng()
                    reac_set = rng.choice(
                        not_used, size=chunk_size, replace=False, shuffle=False
                    )
                    not_used = [_ for _ in not_used if _ not in reac_set]
                else:
                    reac_set = not_used

                # get subset of reactions
                subset_reac = cobra.Model("subset_reac")
                for n in reac_set:
                    subset_reac.add_reactions([universal.reactions[n].copy()])

                # gapfilling
                solution = single_cobra_gapfill(
                    model, subset_reac, medium, namespace, growth_threshold
                )

                if (isinstance(solution, bool) and solution) or (
                    isinstance(solution, list) and len(solution) > 0
                ):
                    break

            if (isinstance(solution, bool) and solution) or (
                isinstance(solution, list) and len(solution) > 0
            ):
                break

    # use the whole reactions content of the universal model at once
    #     not advised for Laptops and PCs with small computing power
    #     may take a long time
    else:
        solution = single_cobra_gapfill(
            model, universal, medium, namespace, growth_threshold
        )

    # if solution is found add new reactions to model
    if isinstance(solution, list) and len(solution) > 0:
        for reac in solution[0]:
            reac.notes["creation"] = "via gapfilling"
        print(
            f"Adding {len(solution[0])} reactions to model to ensure growth on medium {medium.name}."
        )
        model.add_reactions(solution[0])

    return model


def multiple_cobra_gapfill(
    model: cobra.Model,
    universal: cobra.Model,
    media_list: list[Medium],
    namespace: Literal["BiGG"] = "BiGG",
    growth_threshold: float = 0.05,
    iterations: int = 3,
    chunk_size: int = 10000,
) -> cobra.Model:
    """Perform :py:func:`~refinegems.classes.gapfill.cleanup.single_cobra_gapfill` on a list of media.

    Args:
        - model (cobra.Model):
            The model to be gapfilled.
        - universal (cobra.Model):
            The model to use reactions for gapfilling from.
        - media_list (list[Medium]):
            List ofmedia the model is supposed to grow on.
        - growth_threshold (float, optional):
            Growth threshold for the gapfilling.
            Defaults to 0.05.
        - iterations (int, optional):
            Number of iterations for the heuristic version of the gapfilling.
            If 0 or None is given, uses full set of reactions.
            Defaults to 3.
        - chunk_size (int, optional):
            Number of reactions to be used for gapfilling at the same time.
            If None or 0 is given, use full set, not heuristic.
            Defaults to 10000.
    Returns:
        cobra.Model:
            The gapfilled model, if a solution was found.
    """

    for medium in media_list:
        model = cobra_gapfill_wrapper(
            model,
            universal,
            medium,
            namespace,
            iterations,
            chunk_size,
            growth_threshold,
        )

    return model
