# AutoFitLor: Automated Lorentzian Fitting in XSPEC

AutoFitLor is a Tcl script designed to automate the fitting of Lorentzian models to spectral data in XSPEC. It iteratively adds Lorentzian components to fit power density spectra (PDS), real and imaginary parts, phase lags, and coherence spectra, using statistical tests to determine the optimal number of Lorentzians.

## Features

- **Automated Fitting**: Iteratively adds Lorentzians until residuals are random or the model fails statistical tests.
- **Statistical Tests**: Uses the Wald-Wolfowitz runs test and F-test to evaluate model improvements.
- **XSPEC Integration**: Seamlessly integrates with XSPEC for spectral fitting and plotting.
- **Custom Models**: Defines custom Lorentzian models with phase normalization for real and imaginary parts.

## Installation

### Prerequisites

- **XSPEC**: Ensure XSPEC is installed and configured on your system.
- **Tcllib**: The script requires the `math::statistics` package from Tcllib. Install Tcllib and ensure the path is correctly set in the script (default: `/opt/software/tcllib/`).

### Setup

1. Clone the repository:
   ```bash
   git clone https://github.com/[your-username]/autofitlor.git
   cd autofitlor
