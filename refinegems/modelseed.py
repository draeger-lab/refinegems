#!/usr/bin/env python
""" Reports mismatches in charges and formulae based on ModelSEED

Extracts modelseed data from a given tsv file, extracts all metabolites
from a given model. Both lists of metabolites are compared by charge and
formula.
"""
import pandas as pd
import re
import numpy as np


def get_modelseed_compounds(path):
    """extracts compounds from modelseed which have BiGG Ids

    Args:
        path (str): path to modelseed compound definition

    Returns:
        df: table containing modelseed data
    """
    com = pd.read_csv(path, sep='\t')

    def get_bigg_ids(aliases):
        try:
            aliases_list = aliases.split('|')
            bigg = [x[6:] for x in aliases_list if re.search('BiGG: .*', x)]
            return bigg[0]
        except (IndexError, AttributeError):
            return None

    com['BiGG'] = com.apply(lambda row: get_bigg_ids(row['aliases']), axis=1)

    return com.loc[:, ['id', 'name', 'formula', 'mass', 'charge', 'BiGG']].dropna(subset=['BiGG'])


def get_model_charges(model):
    """extracts all metabolites from model

    Args:
        model (cobra-model): model loaded with cobrapy

    Returns:
        df: table containing charges and formulae of model metabolites
    """
    charges = {}
    for metab in model.metabolites:
        charges[metab.id[:-2]] = [metab.charge, metab.formula]

    df_charges = pd.DataFrame.from_dict(charges, orient='index', columns=[
                                        'charge_model', 'formula_model']).reset_index().rename(columns={'index': 'BiGG'})

    return df_charges


def get_modelseed_charges(modelseed_compounds):
    """extract table with BiGG, charges and formulae

    Args:
        modelseed_compounds (df): containing modelseed data

    Returns:
        df: table containing charges and formulae of modelseed metabolites
    """
    modelseed_compounds = modelseed_compounds.loc[:, ['charge', 'BiGG', 'formula']].rename(
        columns={'charge': 'charge_modelseed', 'formula': 'formula_modelseed'})
    lst_col = 'BiGG'
    x = modelseed_compounds.assign(
        **{lst_col: modelseed_compounds[lst_col].str.split(';')})
    df_ms = pd.DataFrame({col: np.repeat(x[col].values, x[lst_col].str.len()) for col in x.columns.difference(
        [lst_col])}).assign(**{lst_col: np.concatenate(x[lst_col].values)})[x.columns.tolist()]
    return df_ms


def compare_model_modelseed(model_charges, modelseed_charges):
    """compares tables with charges / formulae from model & modelseed

    Args:
        model_charges (df): charges and formulae of model metabolites
        modelseed_charges (df): charges and formulae of modelseed metabolites

    Returns:
        df: table containing info whether charges / formulae match
    """
    df_comp = pd.merge(model_charges, modelseed_charges, on='BiGG', how='left')

    def f(x):
        return True if float(x['charge_model']) == x['charge_modelseed'] else False

    def g(x):
        return True if x['formula_model'] == x['formula_modelseed'] else False

    df_comp['charge_match'] = df_comp.apply(f, axis=1)
    df_comp['formula_match'] = df_comp.apply(g, axis=1)

    return df_comp


def get_charge_mismatch(df_comp):
    """extracts metabolites with charge mismatch of model & modelseed

    Args:
        df_comp (df): charge and formula mismatches

    Returns:
        df: table containing metabolites with charge mismatch
    """
    return df_comp.loc[~df_comp['charge_match']].dropna(subset=['charge_modelseed'])


def get_formula_mismatch(df_comp):
    """extracts metabolites with formula mismatch of model & modelseed

    Args:
        df_comp (df): charge and formula mismatches

    Returns:
        df: table containing metabolites with formula mismatch
    """
    return df_comp.loc[~df_comp['formula_match']].dropna(subset=['formula_modelseed'])


def get_compared_formulae(formula_mismatch):
    """compare formula by atom pattern

    Args:
        formula_mismatch (df): table with column containing atom comparison

    Returns:
        df: table containing metabolites with formula mismatch
    """

    def formula_comparison(f1, f2):
        formula_pattern = "[A-Z][a-z]?\d*"
        atom_pattern = '[A-Z][a-z]?'
        atom_number_pattern = '\d+'
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
        lambda row: formula_comparison(row['formula_model'], row['formula_modelseed']), axis=1)

    return formula_mismatch


def modelseed(path, model):
    """Executes all steps to compare model metabolites to modelseed metabolites

    Args:
        path (str): path to modelseed compound definition
        model (cobra-model): model loaded with cobrapy

    Returns:
        tuple: (table with charge mismatches, formula mismatches)
    """
    ms_comp = get_modelseed_compounds(path)
    model_charges = get_model_charges(model)
    modelseed_charges = get_modelseed_charges(ms_comp)
    df_comp = compare_model_modelseed(model_charges, modelseed_charges)
    charge_mismatch = get_charge_mismatch(df_comp)
    formula_mismatch = get_formula_mismatch(df_comp)
    formula_comp = get_compared_formulae(formula_mismatch)
    return charge_mismatch, formula_comp
