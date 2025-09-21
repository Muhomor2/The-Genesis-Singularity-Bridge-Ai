# QadmonAnalyzer: A Universal Data Analysis Tool

**Authors**: [Igor Chechelnitsky]  
**Date**: September 21, 2025  
**Repository**: [https://github.com/Muhomor2/The-Genesis-Singularity-Bridge-Ai/tree/main]  

---

## Overview

`QadmonAnalyzer` is a Python-based tool designed for comprehensive exploratory data analysis (EDA) of structured datasets in CSV format. It provides a suite of analytical methods to uncover patterns in numerical and categorical data, including:

- Detection of significant constants (ontological anchors) in numerical columns.
- Frequency analysis of categorical variables with chi-square testing.
- Time series analysis through autocorrelation and Hurst exponent estimation.
- Identification of repeating sequences in categorical data.

Built with `numpy`, `pandas`, `scipy`, and optional `matplotlib` for visualization, `QadmonAnalyzer` is suitable for researchers and data scientists in fields such as physics, sociology, and bioinformatics. This document describes the tool's functionality, mathematical foundations, usage, and example application.

---

## Features

- **Anchor Detection**: Identifies numerical values close to fundamental constants (e.g., π, e, gravitational acceleration).
- **Categorical Analysis**: Performs frequency analysis and chi-square tests to assess distribution randomness.
- **Time Series Analysis**: Computes autocorrelation and Hurst exponent for numerical series to evaluate trends and persistence.
- **Sequence Detection**: Finds frequent subsequences in categorical data for pattern analysis.
- **Output**: Saves results in JSON format for reproducibility and integration.
- **Visualization**: Optional plotting of distributions and time series using `matplotlib`.

---

## Installation

### Prerequisites
- Python 3.8 or higher
- Required libraries: `numpy`, `pandas`, `scipy`
- Optional: `matplotlib` for visualization

Install dependencies using pip:
```bash
pip install numpy pandas scipy matplotlib
