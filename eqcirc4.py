#!/usr/bin/env python3
"""
Equivalent Circuit Parameter Estimator for Piezoelectric Structures
Enhanced version with logging and configurable plotting options

Author: D. S. Stutts (original), Enhanced by user
Associate Professor of Mechanical Engineering
Missouri University of Science and Technology
Email: stutts@mst.edu

Original release: eqcirc2.py Version 0.1.0 3-29-2015
Enhanced version: 2025

This program calculates equivalent circuit parameters from frequency-impedance
magnitude data stored in the standard HP4294A Impedance Analyzer output format.

Equivalent circuit diagram:
http://web.mst.edu/~stutts/piezoequivcircuit0.png
"""

import numpy as np
from scipy.optimize import leastsq
import matplotlib.pyplot as plt
import matplotlib.patches as mpatch
import logging
import argparse
from pathlib import Path
from datetime import datetime


class PiezoAnalyzer:
    """Analyzes piezoelectric impedance data and fits equivalent circuit model."""

    def __init__(self, log_level=logging.INFO):
        """Initialize analyzer with logging configuration."""
        self.setup_logging(log_level)
        self.logger = logging.getLogger(__name__)

    def setup_logging(self, log_level):
        """Configure logging to file and console."""
        log_filename = f'piezo_analysis_{datetime.now().strftime("%Y%m%d_%H%M%S")}.log'

        logging.basicConfig(
            level=log_level,
            format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
            handlers=[
                logging.FileHandler(log_filename),
                logging.StreamHandler()
            ]
        )

    @staticmethod
    def admittance_model(f, z):
        """
        Calculate admittance from equivalent circuit parameters.

        This is the exact formula from the original eqcirc2.py code.
        Note: 0.2e1 = 2.0, 0.4e1 = 4.0, 0.1e1 = 1.0 in the original

        Parameters:
        -----------
        f : array-like
            Frequency array (Hz)
        z : array-like
            Circuit parameters [C0, R1, L1, C1]
            C0: parallel capacitance (F)
            R1: motional resistance (Ω)
            L1: motional inductance (H)
            C1: motional capacitance (F)

        Returns:
        --------
        Y : array-like
            Admittance (A/V)
        """
        # Exact translation from original:
        # return 0.2e1 * np.pi * f * np.sqrt(0.4e1 * z[0] ** 2*z[3]**2*
        # z[1] ** 2 * np.pi ** 2 * f ** 2 + (
        # -0.4e1 * z[0] * z[3] * z[2] * np.pi ** 2 * f ** 2
        # + z[0]+z[3])**2)*((-0.4e1 * z[3]*z[2]*np.pi ** 2*f**2+0.1e1)**2
        # + 0.4e1*z[1]**2*z[3]**2*np.pi**2*f**2)**(-0.1e1/0.2e1)

        # Broken down for clarity:
        numerator_term1 = 4 * z[0]**2 * z[3]**2 * z[1]**2 * np.pi**2 * f**2
        numerator_term2 = (-4 * z[0] * z[3] * z[2] * np.pi**2 * f**2 + z[0] + z[3])**2
        numerator = numerator_term1 + numerator_term2

        denominator_term1 = (-4 * z[3] * z[2] * np.pi**2 * f**2 + 1)**2
        denominator_term2 = 4 * z[1]**2 * z[3]**2 * np.pi**2 * f**2
        denominator = denominator_term1 + denominator_term2

        return 2 * np.pi * f * np.sqrt(numerator) * (denominator)**(-0.5)

    @staticmethod
    def estimate_C0(Ymin, Ymax, fr, fa):
        """Estimate parallel capacitance."""
        term1 = 2 * (fa**2 - fr**2) * Ymin**2 / (np.pi**2 * fa**4)
        term2 = 2 * np.sqrt((fa**2 - fr**2)**2 * Ymin**4 / (np.pi**4 * fa**8) +
                           4 * Ymin**2 * Ymax**2 / (np.pi**4 * fa**4))
        return np.sqrt(term1 + term2) / 4

    @staticmethod
    def estimate_R1(Ymin, Ymax, fr, fa, C0):
        """Estimate motional resistance."""
        return 1 / np.sqrt(-4 * np.pi**2 * fr**2 * C0**2 + Ymax**2)

    @staticmethod
    def estimate_L1(fr, fa, C0):
        """Estimate motional inductance."""
        return 1 / (4 * np.pi**2 * (fa**2 - fr**2) * C0)

    @staticmethod
    def estimate_C1(fr, fa, C0):
        """Estimate motional capacitance."""
        return (fa**2 / fr**2 - 1) * C0

    def residual(self, z, ydat, f):
        """Calculate residual for least squares optimization."""
        return ydat - self.admittance_model(f, z)

    def load_data(self, filepath):
        """
        Load impedance data from HP4294A format file.

        Parameters:
        -----------
        filepath : str or Path
            Path to data file

        Returns:
        --------
        freq : ndarray
            Frequency array
        impedance : ndarray
            Impedance magnitude array
        """
        self.logger.info(f"Loading data from {filepath}")

        with open(filepath, 'r') as f:
            lines = f.readlines()

        numlines = len(lines)
        nummagpts = (numlines - 1 - 26) // 2

        self.logger.info(f"Total lines: {numlines}, Data points: {nummagpts}")

        freq = []
        impedance = []

        for i, line in enumerate(lines[21:21+nummagpts]):
            try:
                parts = line.split()
                freq.append(float(parts[0]))
                impedance.append(float(parts[1]))
            except (ValueError, IndexError) as e:
                self.logger.warning(f"Skipping line {i+21}: {e}")

        return np.array(freq), np.array(impedance)

    def analyze(self, freq, impedance):
        """
        Perform equivalent circuit parameter estimation.

        Parameters:
        -----------
        freq : ndarray
            Frequency array
        impedance : ndarray
            Impedance magnitude array

        Returns:
        --------
        results : dict
            Dictionary containing all fitted parameters and metrics
        """
        self.logger.info("Starting parameter estimation")

        # Calculate admittance
        admittance = 1 / impedance

        # Find resonance and anti-resonance
        Ymax = np.max(admittance)
        fr = freq[np.argmax(admittance)]
        Ymin = np.min(admittance)
        fa = freq[np.argmin(admittance)]

        self.logger.info(f"Ymax = {Ymax:.4e} at fr = {fr:.4e} Hz")
        self.logger.info(f"Ymin = {Ymin:.4e} at fa = {fa:.4e} Hz")

        # Initial parameter estimates
        C0i = self.estimate_C0(Ymin, Ymax, fr, fa)
        R1i = self.estimate_R1(Ymin, Ymax, fr, fa, C0i)
        L1i = self.estimate_L1(fr, fa, C0i)
        C1i = self.estimate_C1(fr, fa, C0i)

        self.logger.info("Initial estimates:")
        self.logger.info(f"  C0 = {C0i:.4e} F")
        self.logger.info(f"  R1 = {R1i:.4f} Ω")
        self.logger.info(f"  L1 = {L1i:.4e} H")
        self.logger.info(f"  C1 = {C1i:.4e} F")

        # Optimize parameters using Levenberg-Marquardt
        z0 = [C0i, R1i, L1i, C1i]
        output = leastsq(self.residual, z0, args=(admittance, freq), full_output=1)

        C0, R1, L1, C1 = output[0]

        # Calculate derived parameters
        Q = 1 / (R1 * np.sqrt(C1 / L1))
        fr_fit = 1 / (2 * np.pi * np.sqrt(L1 * C1))
        fa_fit = np.sqrt((C0 + C1) / (C0 * C1 * L1)) / (2 * np.pi)

        # Calculate RMS error
        y_model = self.admittance_model(freq, [C0, R1, L1, C1])
        rmserr = np.sqrt(np.sum((admittance - y_model)**2) / len(freq))

        results = {
            'freq': freq,
            'impedance': impedance,
            'admittance': admittance,
            'C0': C0,
            'R1': R1,
            'L1': L1,
            'C1': C1,
            'Q': Q,
            'fr': fr_fit,
            'fa': fa_fit,
            'rmserr': rmserr,
            'y_model': y_model
        }

        self.logger.info("\nOptimized parameters:")
        self.logger.info(f"  fr = {fr_fit:.4e} Hz")
        self.logger.info(f"  fa = {fa_fit:.4e} Hz")
        self.logger.info(f"  C0 = {C0:.4e} F")
        self.logger.info(f"  R1 = {R1:.4f} Ω")
        self.logger.info(f"  L1 = {L1:.4e} H")
        self.logger.info(f"  C1 = {C1:.4e} F")
        self.logger.info(f"  Q = {Q:.2f}")
        self.logger.info(f"  RMS Error = {rmserr:.2e}")

        return results

    def plot_results(self, results, show_model=True, show_annotations=True,
                    journal_style=False, output_file=None, freq_units='Hz'):
        """
        Generate plots of data and fitted model.

        Parameters:
        -----------
        results : dict
            Results dictionary from analyze()
        show_model : bool
            Whether to plot the fitted model
        show_annotations : bool
            Whether to show parameter annotations on plot
        journal_style : bool
            Use journal-ready formatting (larger fonts, no background)
        output_file : str or None
            Output filename for saving plot
        freq_units : str
            Frequency units for plot ('Hz', 'kHz', 'MHz')
        """
        self.logger.info("Generating plot")

        # Determine frequency scaling
        freq_scale = {'Hz': 1.0, 'kHz': 1e-3, 'MHz': 1e-6}
        if freq_units not in freq_scale:
            self.logger.warning(f"Unknown frequency unit '{freq_units}', using 'Hz'")
            freq_units = 'Hz'

        scale = freq_scale[freq_units]
        freq_scaled = results['freq'] * scale

        # Set style based on journal requirements
        if journal_style:
            plt.style.use('seaborn-v0_8-paper')
            fig, ax = plt.subplots(1, figsize=(6, 4), dpi=300)
            fontsize_label = 14
            fontsize_tick = 12
            fontsize_legend = 12
            fontsize_annot = 10
            linewidth = 1.5
            markersize = 2
        else:
            fig, ax = plt.subplots(1, figsize=(10, 6))
            fontsize_label = 12
            fontsize_tick = 12
            fontsize_legend = 12
            fontsize_annot = 12
            linewidth = 2
            markersize = 1

        # Plot model first (behind data) if requested
        if show_model:
            ax.plot(freq_scaled, results['y_model'], 'r-',
                   linewidth=linewidth, label='Fitted Model', zorder=1)

        # Plot data on top
        ax.plot(freq_scaled, results['admittance'], 'ko',
                markersize=markersize, label='Measured Data', zorder=2)

        # Add annotations if requested and not in journal style
        if show_annotations and not journal_style:
            Ymax = np.max(results['admittance'])
            Ymin = np.min(results['admittance'])
            dely = Ymax - Ymin
            delx = freq_scaled[-1] - freq_scaled[0]

            xpos = float(results['fa']) * scale
            noteymax = 0.8 * Ymax
            padx = delx / 100
            legendwidth = delx / 2
            legendheight = dely / 1.75
            linegap = dely / 15

            # Background rectangle
            ax.add_patch(mpatch.Rectangle(
                (xpos - padx, noteymax),
                legendwidth, -legendheight,
                alpha=0.5, facecolor='#ffcccc', zorder=2
            ))

            # Text annotations
            annotations = [
                f"fr = {results['fr']*scale:.4e} {freq_units}",
                f"fa = {results['fa']*scale:.4e} {freq_units}",
                f"C0 = {results['C0']:.4e} F",
                f"R1 = {results['R1']:.4f} Ω",
                f"L1 = {results['L1']:.4e} H",
                f"C1 = {results['C1']:.4e} F",
                f"Q = {results['Q']:.2f}",
                f"RMSErr = {results['rmserr']:.2e}"
            ]

            for i, text in enumerate(annotations):
                ax.annotate(text, xy=(xpos, noteymax - (i+1)*linegap),
                           fontsize=fontsize_annot, zorder=3)

        # Formatting
        ax.set_xlabel(f'Frequency ({freq_units})', fontsize=fontsize_label)
        ax.set_ylabel('Admittance (S)', fontsize=fontsize_label)
        ax.tick_params(labelsize=fontsize_tick)

        # Only use scientific notation for Hz
        if freq_units == 'Hz':
            ax.ticklabel_format(axis='x', style='sci', scilimits=(0, 0))

        ax.grid(True, alpha=0.3)

        legend = ax.legend(loc='upper right', fontsize=fontsize_legend,
                          framealpha=1 if not journal_style else 0.9)

        if not journal_style:
            legend.get_frame().set_facecolor('white')

        plt.tight_layout()

        # Save if requested
        if output_file:
            self.logger.info(f"Saving plot to {output_file}")
            # Determine format from extension, default to SVG for journal style
            if journal_style and not output_file.endswith(('.svg', '.pdf', '.png')):
                output_file = Path(output_file).stem + '.svg'
            plt.savefig(output_file,
                       dpi=300 if journal_style and output_file.endswith('.png') else None,
                       format='svg' if output_file.endswith('.svg') else None,
                       bbox_inches='tight')

        plt.show()


def main():
    """Main execution function with command-line interface."""
    parser = argparse.ArgumentParser(
        description='Piezoelectric Equivalent Circuit Parameter Estimator',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument('input_file', type=str,
                       help='Input data file (HP4294A format)')
    parser.add_argument('--no-model', action='store_true',
                       help='Plot data only without fitted model')
    parser.add_argument('--no-annotations', action='store_true',
                       help='Hide parameter annotations on plot')
    parser.add_argument('--journal', action='store_true',
                       help='Generate journal-ready figure')
    parser.add_argument('--output', '-o', type=str,
                       help='Output plot filename (default: auto-generated)')
    parser.add_argument('--freq-units', type=str, default='Hz',
                       choices=['Hz', 'kHz', 'MHz'],
                       help='Frequency units for plot (default: Hz)')
    parser.add_argument('--log-level', type=str, default='INFO',
                       choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'],
                       help='Logging level (default: INFO)')

    args = parser.parse_args()

    # Initialize analyzer
    analyzer = PiezoAnalyzer(log_level=getattr(logging, args.log_level))

    # Load and analyze data
    try:
        freq, impedance = analyzer.load_data(args.input_file)
        results = analyzer.analyze(freq, impedance)

        # Generate output filename if not specified
        if args.output is None:
            input_path = Path(args.input_file)
            suffix = '_journal' if args.journal else '_model'
            extension = '.svg' if args.journal else '.pdf'
            args.output = input_path.stem + suffix + extension

        # Plot results
        analyzer.plot_results(
            results,
            show_model=not args.no_model,
            show_annotations=not args.no_annotations,
            journal_style=args.journal,
            output_file=args.output,
            freq_units=args.freq_units
        )

        analyzer.logger.info("Analysis completed successfully")

    except Exception as e:
        logging.error(f"Error during analysis: {e}", exc_info=True)
        return 1

    return 0


if __name__ == '__main__':
    exit(main())
