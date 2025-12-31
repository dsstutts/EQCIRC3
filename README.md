# Piezoelectric Equivalent Circuit Parameter Estimator

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.7+](https://img.shields.io/badge/python-3.7+-blue.svg)](https://www.python.org/downloads/)

Enhanced version of the piezoelectric equivalent circuit parameter estimation tool with modern features including logging, configurable plotting, and journal-ready figure generation.

## Overview

This tool analyzes frequency-impedance magnitude data from piezoelectric structures and estimates equivalent circuit parameters using the Levenberg-Marquardt nonlinear least squares algorithm. It processes data from HP4294A Impedance Analyzer output files.

### Equivalent Circuit

The program fits data to the following equivalent circuit model:

![Equivalent Circuit](http://web.mst.edu/~stutts/piezoequivcircuit0.png)

### Calculated Parameters

- **fr** - Resonance frequency (Hz)
- **fa** - Anti-resonance frequency (Hz)
- **C0** - Parallel capacitance (F)
- **R1** - Motional resistance (Ω)
- **L1** - Motional inductance (H)
- **C1** - Motional capacitance (F)
- **Q** - Quality factor (Q = 1/2ζ)
- **RMS Error** - Root mean square fitting error

## Features

✨ **New in Enhanced Version:**
- 🔍 Comprehensive logging system with timestamped log files
- 📊 Toggle options for displaying model overlay
- 📝 Optional parameter annotations on plots
- 📄 Journal-ready figure generation (SVG format, high resolution)
- 🎯 Command-line interface with flexible options
- 🏗️ Clean object-oriented architecture
- ⚠️ Robust error handling
- 🎨 Proper plot layering (data points always visible over fitted model)
- 📏 Configurable frequency units (Hz, kHz, MHz) for plot displays

## Installation

### Requirements

```bash
python >= 3.7
numpy >= 1.19.0
scipy >= 1.5.0
matplotlib >= 3.3.0
```

### Install Dependencies

```bash
pip install numpy scipy matplotlib
```

Or using conda:

```bash
conda install numpy scipy matplotlib
```

## Usage

### Basic Usage

```bash
python eqcirc2_enhanced.py inputfile.txt
```

This will:
- Analyze the impedance data
- Generate a log file with timestamp
- Create a plot with fitted model and parameter annotations
- Save the plot as `inputfile_model.pdf`

### Command-Line Options

```bash
python eqcirc2_enhanced.py [-h] [--no-model] [--no-annotations]
                            [--journal] [--output OUTPUT]
                            [--freq-units {Hz,kHz,MHz}]
                            [--log-level {DEBUG,INFO,WARNING,ERROR}]
                            input_file
```

#### Positional Arguments
- `input_file` - Input data file in HP4294A format (required)

#### Optional Arguments
- `-h, --help` - Show help message and exit
- `--no-model` - Plot measured data only without fitted model overlay
- `--no-annotations` - Hide parameter annotations on plot
- `--journal` - Generate journal-ready figure (SVG, 300 DPI, clean styling)
- `--output OUTPUT, -o OUTPUT` - Specify output plot filename
- `--freq-units {Hz,kHz,MHz}` - Set frequency units for plot axis (default: Hz)
- `--log-level {DEBUG,INFO,WARNING,ERROR}` - Set logging verbosity (default: INFO)

### Examples

#### Standard Analysis with All Features
```bash
python eqcirc2_enhanced.py F5A2.TXT
```
Output: `F5A2_model.pdf` and `piezo_analysis_YYYYMMDD_HHMMSS.log`

#### Data Only (No Model Overlay)
```bash
python eqcirc2_enhanced.py F5A2.TXT --no-model
```
Useful for visualizing raw measurements without the fitted curve.

#### Clean Plot (No Annotations)
```bash
python eqcirc2_enhanced.py F5A2.TXT --no-annotations
```
Shows data and model without the parameter text box.

#### Journal-Ready Figure
```bash
python eqcirc2_enhanced.py F5A2.TXT --journal -o figure1.svg
```
Generates publication-quality SVG with:
- High resolution (300 DPI for raster elements)
- Vector graphics (infinitely scalable)
- Larger, readable fonts
- Clean, professional styling
- No background annotations

#### Custom Output with Debug Logging
```bash
python eqcirc2_enhanced.py F5A2.TXT --log-level DEBUG -o results/analysis.pdf
```
Provides detailed logging information for troubleshooting.

#### Combination of Options
```bash
python eqcirc2_enhanced.py F5A2.TXT --journal --no-annotations -o paper_fig1.svg
```
Creates a minimal, journal-ready figure with just data and model curves.

#### Frequency Units
```bash
# Display frequency in kHz (useful for 10-100 kHz range)
python eqcirc2_enhanced.py F5A2.TXT --freq-units kHz

# Display frequency in MHz (useful for MHz range measurements)
python eqcirc2_enhanced.py F5A2.TXT --freq-units MHz --journal -o high_freq.svg
```
Changes the frequency axis units and updates all displayed values accordingly. Internal calculations remain in Hz for accuracy.

## Input File Format

The program expects data files in HP4294A Impedance Analyzer format:
- Header information in first 21 lines
- Data section with frequency and impedance magnitude columns
- Standard format: `frequency(Hz)  impedance(Ω)`

Example:
```
[Header lines 1-21]
5.000000e+04  1.234567e+03
5.100000e+04  1.245678e+03
...
```

## Output Files

### Log File
- Filename: `piezo_analysis_YYYYMMDD_HHMMSS.log`
- Contains: Analysis parameters, fitted values, and processing information
- Levels: INFO (default), DEBUG, WARNING, ERROR

### Plot File
- Default naming:
  - Standard mode: `inputfile_model.pdf`
  - Journal mode: `inputfile_journal.svg`
- Custom naming: Use `-o` or `--output` flag
- Supported formats: SVG, PDF, PNG (auto-detected from extension)
- Plot layering: Measured data points are always drawn on top of the fitted model curve for maximum visibility

### Console Output
All fitted parameters are printed to console:
```
fr =  5.4321e+04 Hz
fa =  5.5678e+04 Hz
C0 =  1.234e-09 F
R1 =  12.3456 Ω
L1 =  5.678e-02 H
C1 =  9.876e-12 F
Q =  123.45
RMS Error = 1.23e-05
```

**Note**: Frequency values in console output are always in Hz regardless of the `--freq-units` setting. The `--freq-units` flag only affects the plot display for better visualization.

## Algorithm Details

The program uses a two-stage approach:

1. **Initial Estimation**: Analytical formulas provide approximate circuit parameters based on resonance and anti-resonance frequencies
2. **Optimization**: Levenberg-Marquardt algorithm refines parameters using nonlinear least squares fitting

The admittance model used for fitting is derived from the equivalent circuit's impedance transfer function.

## Publication Guidelines

### Citing This Software

If you use this software in your research, please cite:

```bibtex
@Misc{eqcirc2_2015,
  author = {Stutts, D. S.},
  title = {{eqcirc2.py}: {Equivalent Circuit Parameter Estimator for Piezoelectric Structures}},
  howpublished = {\url{https://github.com/MSTESG/EQCIRC1.git}},
  year = {2015}
}
```

### Journal Figure Best Practices

When preparing figures for publication:

1. Use `--journal` flag for optimized formatting
2. Save as SVG for vector quality: `-o figure.svg`
3. Use `--no-annotations` for cleaner appearance
4. Select appropriate frequency units with `--freq-units`:
   - Use `kHz` for typical piezoelectric resonators (10-100 kHz range)
   - Use `MHz` for high-frequency devices (>1 MHz)
   - Use `Hz` with scientific notation for mixed ranges
5. Plots automatically layer data points over fitted curves for clarity
6. Consider journal's specific requirements for:
   - Font sizes (currently: 14pt labels, 12pt ticks, 10pt annotations)
   - Figure dimensions (currently: 6×4 inches)
   - Color schemes (currently: black data points, red fitted model)

## Troubleshooting

### Common Issues

**Issue**: `FileNotFoundError: [Errno 2] No such file or directory`
- **Solution**: Check that input file path is correct and file exists

**Issue**: Plot appears but doesn't save
- **Solution**: Ensure write permissions in output directory

**Issue**: `ValueError: could not convert string to float`
- **Solution**: Verify input file format matches HP4294A specification

**Issue**: Poor fit quality (high RMS error)
- **Solution**:
  - Check data quality and frequency range
  - Ensure resonance/anti-resonance peaks are clearly defined
  - Verify data covers appropriate frequency span around resonance

### Debug Mode

For detailed troubleshooting information:
```bash
python eqcirc2_enhanced.py inputfile.txt --log-level DEBUG
```

This provides:
- Data loading details
- Initial parameter estimates
- Optimization progress
- Detailed error messages

## Development

### Project Structure

```
eqcirc2_enhanced.py
├── PiezoAnalyzer (class)
│   ├── __init__()
│   ├── setup_logging()
│   ├── admittance_model()
│   ├── estimate_C0/R1/L1/C1()
│   ├── load_data()
│   ├── analyze()
│   └── plot_results()
└── main()
```

### Contributing

Contributions are welcome! Areas for enhancement:
- Additional equivalent circuit models
- Support for other impedance analyzer formats
- Multi-resonance fitting
- Automated quality assessment
- GUI interface
- Additional frequency unit options (GHz, etc.)

## License

This code is released under the MIT License.

Original code copyright (c) 2015 D. S. Stutts
Enhanced version copyright (c) 2025

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

## Contact

**Original Author:**
D. S. Stutts
Associate Professor of Mechanical Engineering
Missouri University of Science and Technology
Email: stutts@mst.edu

## Acknowledgments

- Original eqcirc2.py developed by Dr. D. S. Stutts at Missouri S&T
- Enhanced version with modern Python features and publication tools
- Built using NumPy, SciPy, and Matplotlib

---

**Version**: 2.0 Enhanced
**Last Updated**: December 2025
**Python Compatibility**: 3.7+
