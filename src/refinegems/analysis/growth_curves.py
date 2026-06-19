#!/usr/bin/env python
"""Provides functions to parse growth plate data, normalize replicates, and calculate growth curves using the Gompertz model."""

__author__ = "Famke Baeuerle, Gwendolyn O. Döbel, Carolin Brune, cb-Hades"

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import warnings
from pathlib import Path
from typing import Union, Dict, List, Tuple

def gompertz_model(t: np.ndarray, A: float, mu_m: float, lag: float) -> np.ndarray:
    """
    Modified Gompertz model for microbial growth (Zwietering et al., 1990).
    
    Args:
        - t (np.ndarray): Time array
        - A (float): Maximum growth (Upper Asymptote)
        - mu_m (float): Maximum absolute growth rate
        - lag (float): Lag time
        
    Returns:
        np.ndarray: Calculated OD values for the given time points.
    """
    return A * np.exp(-np.exp((mu_m * np.e / A) * (lag - t) + 1.0))

def fit_growth_curve(time_data: np.ndarray, od_data: np.ndarray) -> Tuple[list, float]:
    """
    Fits the Gompertz model to a growth curve using non-linear least squares.
    
    Args:
        - time_data (np.ndarray): Array of time points.
        - od_data (np.ndarray): Array of measured Optical Density values.
        
    Returns:
        Tuple[list, float]: The optimized parameters [A, mu_m, lag] and the R-squared value.
    """
    # Sanitize data: replace negative ODs (plate reader artifacts) with 0
    od_data = np.clip(od_data, 0, None)
    
    # Smart initial guesses
    guess_A = np.max(od_data)
    guess_mu = np.max(np.gradient(od_data, time_data)) if len(time_data) > 1 else 0.1
    guess_lag = time_data[np.argmax(np.gradient(od_data, time_data))] / 2 if len(time_data) > 1 else 1.0
    
    p0 = [guess_A, guess_mu, guess_lag]
    
    try:
        popt, _ = curve_fit(gompertz_model, time_data, od_data, p0=p0, bounds=(0, np.inf), maxfev=10000)
        
        # Calculate R-squared for Goodness of Fit
        residuals = od_data - gompertz_model(time_data, *popt)
        ss_res = np.sum(residuals**2)
        ss_tot = np.sum((od_data - np.mean(od_data))**2)
        r_squared = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0.0
        
        return list(popt), r_squared
        
    except Exception as e:
        warnings.warn(f"Curve fitting failed: {e}")
        return [np.nan, np.nan, np.nan], np.nan

def process_plate_data(file_path: Union[str, Path], time_col: str = "Time", condition_mapping: Dict[str, List[str]] = None) -> pd.DataFrame:
    """
    Reads plate reader data, normalizes replicates if provided, and calculates growth parameters.
    
    Args:
        - file_path (Union[str, Path]): Path to the CSV or Excel file containing the data.
        - time_col (str): Name of the column containing time points.
        - condition_mapping (Dict[str, List[str]], optional): Dictionary mapping a condition name 
          to a list of well column names (e.g., {"WT_Glucose": ["A1", "A2", "A3"]}). 
          If None, fits every well individually.
        
    Returns:
        pd.DataFrame: Table of computed parameters (A, mu_m, lag, R_squared) for each condition/well.
    """
    file_path = Path(file_path)
    
    # Handle different input formats
    if file_path.suffix.lower() == '.csv':
        df = pd.read_csv(file_path)
    elif file_path.suffix.lower() in ['.xls', '.xlsx']:
        df = pd.read_excel(file_path)
    else:
        raise ValueError("Unsupported file format. Please provide a .csv or .xlsx file.")
        
    time_data = df[time_col].values
    results = []
    
    # If no mapping is provided, treat every column (except time) as a separate condition
    if condition_mapping is None:
        wells = [col for col in df.columns if col != time_col]
        condition_mapping = {well: [well] for well in wells}

    for condition, wells in condition_mapping.items():
        # Normalise biological/technical replicates by averaging them
        replicate_data = df[wells].values
        od_mean = np.mean(replicate_data, axis=1)
        od_std = np.std(replicate_data, axis=1) if len(wells) > 1 else np.zeros_like(od_mean)
        
        popt, r_sq = fit_growth_curve(time_data, od_mean)
        
        results.append({
            "Condition": condition,
            "Replicates_Averaged": len(wells),
            "Asymptote (A)": popt[0],
            "Max Growth Rate (mu_m)": popt[1],
            "Lag Time (lag)": popt[2],
            "R_squared": r_sq
        })
        
    return pd.DataFrame(results)

def plot_growth_curve(time_data: np.ndarray, od_data: np.ndarray, popt: list, condition_name: str = "Condition", od_std: np.ndarray = None) -> matplotlib.figure.Figure:
    """
    Visualizes the experimental data (with optional error bars) against the fitted Gompertz curve.
    """
    fig = plt.figure(figsize=(8, 5))
    ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
    
    if od_std is not None and np.any(od_std > 0):
        ax.errorbar(time_data, od_data, yerr=od_std, fmt='o', capsize=3, alpha=0.6, label=f'{condition_name} (Mean $\\pm$ SD)')
    else:
        ax.plot(time_data, od_data, 'o', alpha=0.6, label=f'{condition_name} (Raw Data)')
    
    if not np.isnan(popt[0]):
        t_fit = np.linspace(min(time_data), max(time_data), 100)
        y_fit = gompertz_model(t_fit, *popt)
        ax.plot(t_fit, y_fit, '-', color='red', linewidth=2, label=f'Gompertz Fit\n$A$={popt[0]:.2f}, $\\mu_m$={popt[1]:.2f}, lag={popt[2]:.2f}')
        
    ax.set_xlabel('Time')
    ax.set_ylabel('Optical Density (OD)')
    ax.set_title(f'Growth Curve Fit: {condition_name}')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    return fig