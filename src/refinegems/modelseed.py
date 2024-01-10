#!/usr/bin/env python
""" Reports mismatches in charges and formulae based on ModelSEED

Extracts ModelSEED data from a given tsv file, extracts all metabolites from a given model. Both lists of metabolites are compared by charge and formula.
"""
import pandas as pd
import re
import numpy as np
from refinegems.io import load_a_table_from_database
from cobra import Model as cobraModel

__author__ = "Famke Baeuerle and Jan-Philipp Leusch"


def get_modelseed_compounds() -> pd.DataFrame:
    """Extracts compounds from ModelSEED which have BiGG Ids

    Returns:
        pd.DataFrame: Table containing ModelSEED data
    """
    # Get only rows where BiGG is contained
    com = load_a_table_from_database("SELECT *, INSTR(aliases, 'BiGG:') bigg FROM modelseed_compounds WHERE bigg > 0")

    def get_bigg_ids(aliases):
        try:
            aliases_list = aliases.split('|')
            bigg = [x[6:] for x in aliases_list if re.search('BiGG: .*', x)]
            return bigg[0]
        except (IndexError, AttributeError):
            return None

    com['BiGG'] = com.apply(lambda row: get_bigg_ids(row['aliases']), axis=1)

    return com.loc[:, ['id', 'name', 'formula', 'mass', 'charge', 'BiGG']]


def get_model_charges(model: cobraModel) -> pd.DataFrame:
    """Extracts all metabolites from model

    Args:
        - model (cobraModel): Model loaded with COBRApy

    Returns:
        pd.DataFrame: Table containing charges and formulae of model metabolites
    """
    charges = {}
    for metab in model.metabolites:
        charges[metab.id[:-2]] = [metab.charge, metab.formula]

    df_charges = pd.DataFrame.from_dict(
        charges,
        orient='index',
        columns=[
            'charge_model',
            'formula_model']).reset_index().rename(
        columns={
            'index': 'BiGG'})

    return df_charges


def get_modelseed_charges(modelseed_compounds: pd.DataFrame) -> pd.DataFrame:
    """Extract table with BiGG, charges and formulae

    Args:
        - modelseed_compounds (pd.DataFrame): ModelSEED data. Output from get_modelseed_compounds.

    Returns:
        pd.DataFrame: Table containing charges and formulae of ModelSEED metabolites
    """
    modelseed_compounds = modelseed_compounds.loc[:, ['charge', 'BiGG', 'formula']].rename(
        columns={'charge': 'charge_modelseed', 'formula': 'formula_modelseed'})
    lst_col = 'BiGG'
    x = modelseed_compounds.assign(
        **{lst_col: modelseed_compounds[lst_col].str.split(';')})
    df_ms = pd.DataFrame({col: np.repeat(x[col].values, x[lst_col].str.len()) for col in x.columns.difference(
        [lst_col])}).assign(**{lst_col: np.concatenate(x[lst_col].values)})[x.columns.tolist()]
    return df_ms


def compare_model_modelseed(model_charges: pd.DataFrame, modelseed_charges: pd.DataFrame) -> pd.DataFrame:
    """Compares tables with charges / formulae from model & modelseed

    Args:
        - model_charges (pd.DataFrame): Charges and formulae of model metabolites. Output of get_model_charges.
        - modelseed_charges (pd.DataFrame): Charges and formulae of ModelSEED metabolites. Output of get_modelseed_charges.

    Returns:
        pd.DataFrame: Table containing info whether charges / formulae match
    """
    df_comp = pd.merge(model_charges, modelseed_charges, on='BiGG', how='left')

    def f(x):
        return True if float(
            x['charge_model']) == x['charge_modelseed'] else False

    def g(x):
        return True if x['formula_model'] == x['formula_modelseed'] else False

    df_comp['charge_match'] = df_comp.apply(f, axis=1)
    df_comp['formula_match'] = df_comp.apply(g, axis=1)

    return df_comp


def get_charge_mismatch(df_comp: pd.DataFrame) -> pd.DataFrame:
    """Extracts metabolites with charge mismatch of model & modelseed

    Args:
        df_comp (pd.DataFrame): Charge and formula mismatches. Output from compare_model_modelseed.

    Returns:
        pd.DataFrame: Table containing metabolites with charge mismatch
    """
    return df_comp.loc[~df_comp['charge_match']].dropna(
        subset=['charge_modelseed'])


def get_formula_mismatch(df_comp: pd.DataFrame) -> pd.DataFrame:
    """Extracts metabolites with formula mismatch of model & modelseed

    Args:
        df_comp (pd.DataFrame): Charge and formula mismatches. Output from compare_model_modelseed.

    Returns:
        pd.DataFrame: Table containing metabolites with formula mismatch
    """
    return df_comp.loc[~df_comp['formula_match']].dropna(
        subset=['formula_modelseed'])


def get_compared_formulae(formula_mismatch: pd.DataFrame) -> pd.DataFrame:
    """Compare formula by atom pattern

    Args:
        formula_mismatch (pd.DataFrame): Table with column containing atom comparison. Output from get_formula_mismatch.

    Returns:
        pd.DataFrame: table containing metabolites with formula mismatch
    """

    def formula_comparison(f1, f2):
        # from Jan Leusch
        formula_pattern = "[A-Z][a-z]?\\d*"
        atom_pattern = '[A-Z][a-z]?'
        atom_number_pattern = '\\d+'
        difference = {}
        f1_dict = {}
        f2_dict = {}
        f1 = re.findall(formula_pattern, f1)
        f2 = re.findall(formula_pattern, f2)
        for p in f1:
            key = re.findall(atom_pattern, p)[0]
            value = re.findall(atom_number_pattern, p)
            if not value:
                value = 1
            else:
                value = (int)(value[0])
            f1_dict[key] = value
        for q in f2:
            key = re.findall(atom_pattern, q)[0]
            value = re.findall(atom_number_pattern, q)
            if not value:
                value = 1
            else:
                value = (int)(value[0])
            f2_dict[key] = value
        difference = f1_dict
        if not f2_dict:
            return difference
        else:
            for q in f2_dict:
                if q in difference:
                    difference[q] -= f2_dict[q]
                else:
                    difference[q] = -f2_dict[q]
        for a in list(difference):
            if difference[a] == 0:
                difference.pop(a)
        return difference

    formula_mismatch['formula_comparison'] = formula_mismatch.apply(
        lambda row: formula_comparison(
            row['formula_model'], row['formula_modelseed']), axis=1)

    return formula_mismatch


def compare_to_modelseed(model: cobraModel) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Executes all steps to compare model metabolites to ModelSEED metabolites

    Args:
        - model (cobraModel): Model loaded with COBRApy

    Returns:
        tuple: Tables with charge (1) & formula (2) mismatches
            (1) pd.DataFrame: Table with charge mismatches 
            (2) pd.DataFrame: Table with formula mismatches
    """
    ms_comp = get_modelseed_compounds()
    model_charges = get_model_charges(model)
    modelseed_charges = get_modelseed_charges(ms_comp)
    df_comp = compare_model_modelseed(model_charges, modelseed_charges)
    charge_mismatch = get_charge_mismatch(df_comp)
    formula_mismatch = get_formula_mismatch(df_comp)
    formula_comp = get_compared_formulae(formula_mismatch)
    return charge_mismatch, formula_comp
