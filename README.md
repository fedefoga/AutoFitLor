# AutoFitLor: Automated Lorentzian Fitting in XSPEC

AutoFitLor is a Tcl script designed to automate the fitting of Lorentzian models in XSPEC. It iteratively adds Lorentzian components to fit power density spectra (PDS), real and imaginary parts, phase lags, and coherence spectra, using statistical tests to determine the optimal number of Lorentzians. Refer to Méndez et. al (https://academic.oup.com/mnras/article/527/3/9405/7476004) for details on the model.


## Features

- **Automated Fitting**: Iteratively adds Lorentzians until residuals are random or the model fails statistical tests.
- **Statistical Tests**: Uses the Wald-Wolfowitz runs test and F-test to evaluate model improvements.
- **XSPEC Integration**: Seamlessly integrates with XSPEC for spectral fitting and plotting.

## Installation
 - No installation required.

### Prerequisites

- **XSPEC**: Ensure XSPEC is installed and configured on your system.
- **Tcllib**: The script requires the `math::statistics` package from Tcllib. Install Tcllib and ensure the path is correctly set in the script.

### Setup

1. Clone the repository.
2. Start XSPEC and load the following data in order.
      1. Powerspectra in a soft band (eg 0.3-2 keV)
      2. Powerspectra in a hard band (eg 2-12 keV)
      3. Real part of the Crossspectrum
      4. Imaginary part of the crossspectrum
      5. Phase lag spectrum
      6. Intrinsic coherence function
3. Load AutoFitLor (`@AutoFitLor.tcl`).
4. Run the process: `autofitlor`.
