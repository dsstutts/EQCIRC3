####################--
# File: eqcirc4.py
# Equivalent Circuit Parameter Estimator for Piezoelectric Structures
# Author: D. S. Stutts
# Associate Professor of Mechanical Engineering
# 282 Toomey Hall
# 400 W. 13th St.
# Rolla, MO 65409-0050
# Email: stutts@mst.edu
# Original release: eqcirc1.py Version 0.1.0 3-29-2015
# Modified and renamed to eqcirc2.py to input impedance 4-17-2015:DSS
# Now also automatically determines the number of points to process.
#
# Modified for partial Python 3 compatibility, moved results printing
# into functions, refined plot formatting, and added
# flag to set the saved figure format. 4-26-2016: DSS
# Added time stamp on figures. 7-7-2016 Dalton Stover
# Added phase model. 7-8-2016: DSS
# Moved time stamp plot plot title and appended it to the plotfile, and
# added LaTeX formatting to the figure annotations 7-13-2016: DSS
# Added error bar functionality, and changed rmserr to calculate the
# more conservative sqrt(variance/(n-4)) estimate --
# reflecting the four degrees of freedom in the model.
# The error bars are ±rmserr centered on the model
# evaluated at the data point frequency. 1-2-2017: DSS
# Vers 4.0 Enhanced for complex impedance data support 
# and phase data plotting: 6-9-2025 DSS
####################--

"""
    
    This program calculates the equivalent
    circuit parameters from frequency-impedance magnitude
    magnitude data stored in the standard
    HP4294A Impedance Analyser output data
    format.
    
    An image of the equivalent circuit may be
    found here:
    http://web.mst.edu/~stutts/piezoequivcircuit0.png
    
    The program first calculates the approximate
    equivalent circuit parameters for a single
    resonance-antiresonance frequency pair.
    It then uses the Levenberg-Marquardt (LM)
    nonlinear least squares algorithm to
    optimize the model in the least squares
    sense.  The LM algorithm is invoked
    via a call to leastsq(rez, z0, args=(yy, xx),
    full_output=1) from the scipy.optimize library.
    
    See: http://docs.scipy.org/doc/scipy-0.14.0/
    reference/generated/scipy.optimize.leastsq.html
    for more information.
    
    eqcirc3.py calculates the following outputs stdout:
    
    (1) fr (the resonance frequency)
    (2) fa (the anti-resonance frequency)
    (3) C0 (the parallel capacitance)
    (4) R1 (the motional resistance)
    (5) L1 (the motional inductance)
    (6) C1 (the motional capacitance)
    (7) Q (the series R1L1C1 resonance
    quality factor = 1/2zeta)
    (8) RMS Deviation
    
    A graph of the data and model is also produced.
    
    Example call: python eqcirc3.py inputdatafile.txt
    
    The graph may be saved in PNG format, and the text
    may be redirected from stdout to a file like so:
    
    python eqcirc3.py inputdatafile.txt > outdata.txt
    
    Enhanced features in this version:
    - Supports complex impedance data (real + imaginary components)
    - Automatically computes magnitude from complex data
    - Reads and plots both magnitude and phase data
    - Improved data parsing for multiple trace formats
    
    # This code is copyrighted by the author, but released under the MIT
    # license:
    
    Copyright (c) 2015 -- eqcirc3.py
    
    S&T and the University of Missouri Board of Curators
    license to you the right to use, modify, copy, and distribute this
    code subject to the MIT license:
    
    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:
    
    The above copyright notice and this permission notice shall be included
    in all copies or substantial portions of the Software.
    
    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
    THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
    FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
    DEALINGS IN THE SOFTWARE.
    
    The author kindly requests that any publications benefiting from the use
    of this software include the following citation:
    
    @Misc{eqcirc3_2015,
    author =   {Stutts, D. S.},
    title = {{eqcirc3.py}: {Equivalent Circuit Parameter Estimator
    for Piezoelectric Structures.}},
    howpublished = {\\url{https://github.com/MSTESG/EQCIRC3.git}},
    year = {2016}}
    
    """

import sys
from scipy.optimize import leastsq
import numpy as np
from numpy import array, sqrt
import matplotlib.pyplot as plt
import time

# Test for Python version:
cur_version = sys.version_info

# Initialize some lists:
ydat = []
x = []
xx = []
yy = []
zdat = []
phase_data = []
phase_freq = []
f = 0.0

# Define functions:
def y(f, z):  # Admittance model
    return 0.2e1 * np.pi * f * np.sqrt(0.4e1 * z[0] ** 2*z[3]**2*
    2 * np.pi ** 2 * f ** 2 + (-0.4e1 *z[0]*z[3]*z[2] * np.pi**2*f**2
    + z[0]+z[3])**2)*((-0.4e1 * z[3]*z[2]*np.pi ** 2*f**2+0.1e1)**2
    + 0.4e1*z[1]**2*z[3]**2*np.pi**2*f**2)**(-0.1e1/0.2e1)

# Real admittance:
Y_R = lambda f,cg,cg5,cg3,cg1: cg5 * cg1 ** 2 / (16 * cg ** 2 *
cg1 ** 2 * cg3 ** 2 * np.pi ** 4 * f ** 4 + 4 * cg ** 2 * cg1 ** 2 *
np.pi ** 2 * cg5 ** 2 * f ** 2 - 8 * cg ** 2 * cg1 * cg3 * np.pi ** 2
* f ** 2 - 8 * cg * cg1 ** 2 * cg3 * np.pi ** 2 * f ** 2 + cg ** 2 +
2 * cg * cg1 + cg1 ** 2)

# Imaginary admittance:
Y_I = lambda f,cg,cg5,cg3,cg1: -0.1e1 / np.pi / f * (16 * cg * cg1 ** 2
* cg3 ** 2 * np.pi ** 4 * f ** 4 + 4 * cg * cg1 ** 2 * np.pi ** 2
* cg5 ** 2 * f ** 2 - 8 * cg * cg1 * cg3 * np.pi ** 2 * f ** 2 -
4 * cg1 ** 2 * cg3 * np.pi ** 2 * f ** 2 + cg + cg1) / (16 * cg ** 2 *
cg1 ** 2 * cg3 ** 2 * np.pi ** 4 * f ** 4 + 4 * cg ** 2 * cg1 ** 2 *
np.pi ** 2 * cg5 ** 2 * f ** 2 - 8 * cg ** 2 * cg1 * cg3 *
np.pi ** 2 * f ** 2 - 8 * cg * cg1 ** 2 * cg3 * np.pi ** 2 * f ** 2
+ cg ** 2 + 2 * cg * cg1 + cg1 ** 2) / 2

phi = lambda f, cg,cg5,cg3,cg1:180*np.arctan2(Y_I(f,cg,cg5,cg3,cg1),
Y_R(f,cg,cg5,cg3,cg1))/np.pi

def C0_i(Ymin, Ymax, fr, fa):  # Parallel capacitance estimate
    return np.sqrt(0.2e1*(fa ** 2 - fr**2)*Ymin**2/np.pi**2/fa**4
    + 0.2e1*np.sqrt((fa ** 2 - fr ** 2)**2/np.pi**4/fa**8*Ymin**4
    + 0.4e1*Ymin**2*Ymax**2/np.pi**4/fa**4))/0.4e1

def R1_i(Ymin, Ymax, fr, fa, C0):  # Motional resistance estimate
    return (-0.4e1*np.pi**2*fr**2*C0**2+Ymax**2)**(-0.1e1/0.2e1)

def L1_i(fr, fa, C0):  # Motional inductance estimate
    return 0.1e1 / np.pi ** 2 / (fa ** 2 - fr ** 2) / C0 / 0.4e1

def C1_i(fr, fa, C0):  # Motional capacitance estimate
    return (fa ** 2 / fr ** 2 - 1) * C0

def rez(z, ydat, f):  # Residual function
    return ydat - y(f, z)

def parse_data_line(line):
    """Parse a data line and return frequency, real, and imaginary parts"""
    try:
        parts = line.strip().split('\t')
        if len(parts) >= 3:
            freq = float(parts[0])
            real_part = float(parts[1])
            imag_part = float(parts[2])
            return freq, real_part, imag_part
        elif len(parts) >= 2:
            freq = float(parts[0])
            real_part = float(parts[1])
            return freq, real_part, 0.0
    except (ValueError, IndexError):
        pass
    return None

def is_data_line(line):
    """Check if a line contains numerical data"""
    line = line.strip()
    if not line or line.startswith('"') or line.startswith('Frequency'):
        return False
    try:
        parts = line.split('\t')
        float(parts[0])  # Try to parse first element as frequency
        return True
    except (ValueError, IndexError):
        return False

def py3print():
    '''Print results in Python 3.x format'''
    print("Ymax = ", Ymax, " at fr = ", fr, "\n")
    print("Ymin = ", Ymin, " at fa = ", fa, "\n")
    print("fr = ", fr, "\n")
    print("fa = ", fa, "\n")
    # Initial estimates:
    print("C0i = ", C0i,"\n")
    print("R1i = ", R1i,"\n")
    print("L1i = ", L1i,"\n")
    print("C1i = ", C1i,"\n")
    print("Qi = ", Qi,"\n")
    # Optimal estimates:
    print("C0 = ", C0, "\n")
    print("R1 = ", R1, "\n")
    print("L1 = ", L1, "\n")
    print("C1 = ", C1, "\n")
    print("Q = ", Q, "\n")
    print("k31 = ", k31, "\n")
    print("RMS Deviation = ", rmserr,"\n")

# Set the desired resolution:
res = 300  # Use a larger value for PNG
plottype = 'PNG'  # or 'EPS'

# Input data file on command line:
if len(sys.argv) < 2:
    print("Usage: python eqcirc3_improved.py inputdatafile.txt")
    sys.exit(1)

infile = sys.argv[1]

# Read and parse the data file
try:
    with open(infile, "r") as data:
        lines = data.readlines()
except FileNotFoundError:
    print(f"Error: Could not find file {infile}")
    sys.exit(1)

# Parse impedance magnitude data (TRACE A)
impedance_data = []
phase_data_raw = []
current_trace = None
reading_data = False

for i, line in enumerate(lines):
    line = line.strip()
    
    # Detect trace sections
    if 'TRACE: A' in line:
        current_trace = 'A'
        reading_data = False
        continue
    elif 'TRACE: B' in line:
        current_trace = 'B'
        reading_data = False
        continue
    
    # Skip header lines and start reading data after frequency header
    if 'Frequency' in line and 'Data Trace' in line:
        reading_data = True
        continue
    
    # Read numerical data
    if reading_data and is_data_line(line):
        parsed = parse_data_line(line)
        if parsed:
            freq, real_part, imag_part = parsed
            if current_trace == 'A':
                # For impedance data, compute magnitude from complex components
                magnitude = np.sqrt(real_part**2 + imag_part**2)
                impedance_data.append((freq, magnitude))
            elif current_trace == 'B':
                # For phase data, use real part (phase is typically in real component)
                phase_data_raw.append((freq, real_part))

# Convert to arrays
if impedance_data:
    x = [item[0] for item in impedance_data]
    zdat = [item[1] for item in impedance_data]
    xx = array(x)
    zin = array(zdat)
    yy = 1/zin
else:
    print("Error: No impedance data found in file")
    sys.exit(1)

if phase_data_raw:
    phase_freq = array([item[0] for item in phase_data_raw])
    phase_data = array([item[1] for item in phase_data_raw])
    print(f"Found {len(phase_data)} phase data points")
else:
    print("Warning: No phase data found in file")

# Locate Ymax, Ymin, and initial guesses for fr and fa:
zmin = min(zdat)
zminidx = zdat.index(zmin)
zmax = max(zdat)
zmaxidx = zdat.index(zmax)
Ymax = max(yy)
Ymin = min(yy)
fr = x[zminidx]
fa = x[zmaxidx]

# Estimate initial parameter values:
C0i = C0_i(Ymin, Ymax, fr, fa)
R1i = R1_i(Ymin, Ymax, fr, fa, C0_i(Ymin, Ymax, fr, fa))
L1i = L1_i(fr, fa, C0_i(Ymin, Ymax, fr, fa))
C1i = C1_i(fr, fa, C0_i(Ymin, Ymax, fr, fa))
Qi = 1/(R1i*np.sqrt(C1i/L1i))

# Create initial guess array:
z0 = [C0i, R1i, L1i, C1i]

# Find the best values:
output = leastsq(rez, z0, args=(yy, xx), full_output=1)

C0 = output[0][0]
R1 = output[0][1]
L1 = output[0][2]
C1 = output[0][3]
Q = 1 / (R1 * np.sqrt(C1 / L1))
fr = 1 / np.sqrt(L1 * C1) / 0.2e1 / np.pi
fa = np.sqrt((C0 + C1) / C0 / C1 / L1) / np.pi / 0.2e1

# Put the optimal values in a list:
coeffs = [C0, R1, L1, C1]

# Calculate RMS error:
var = np.inner(output[2]['fvec'],output[2]['fvec'])
rmserr = sqrt(var/(len(output[2]['fvec'])-4))
k31 = np.sqrt((fa**2-fr**2)/fa**2)

# Print the results to stdout:
py3print()

# Calculate plot annotation positions:
delx = (fa-fr)/10.0
dely = (Ymax-Ymin)/20.0
noteymax = 0.65*Ymax

# Add date and time:
loctime = time.asctime(time.localtime(time.time()))
print("Date and Time =", loctime, "\n")

# Create the plot
plt.figure(figsize=(8,7), dpi=res)

# Top subplot: Admittance magnitude
plt.subplot(211)
plt.plot(xx, y(xx, coeffs), 'r-', label='model', linewidth=2)
plt.plot(xx, yy, 'go', label='data', markersize=4)
plt.errorbar(xx, y(xx, coeffs), yerr=rmserr, linestyle='None', alpha=0.5)

# Annotations
plt.annotate(r"$f_r$ = "+'{: 3.3e}'.format(fr),xy=(fa-delx,noteymax))
plt.annotate(r"$f_a$ = "+'{: 3.3e}'.format(fa),xy=(fa-delx,noteymax-1.4*dely))
plt.annotate(r"$C_0$ = "+'{: 3.3e}'.format(C0),xy=(fa-delx,noteymax-2.8*dely))
plt.annotate(r"$R_1$ = "+'{: 3.3e}'.format(R1),xy=(fa-delx,noteymax-4.4*dely))
plt.annotate(r"$L_1$ = "+'{: 3.3e}'.format(L1),xy=(fa-delx,noteymax-5.8*dely))
plt.annotate(r"$C_1$ = "+'{: 3.3e}'.format(C1),xy=(fa-delx,noteymax-7.2*dely))
plt.annotate('Q = '+'{: 3.3e}'.format(Q),xy=(fa-delx,noteymax-8.6*dely))
plt.annotate('RMS Dev. = '+'{: 3.2e}'.format(rmserr),xy=(fa-delx,noteymax-10*dely))
plt.annotate('$k_{31}$ = '+'{: 3.2e}'.format(k31),xy=(fa-delx,noteymax-11.4*dely))

plt.suptitle('Data File = '+infile+':  '+loctime)
legend = plt.legend(loc='upper right', shadow=True, fontsize='large')
plt.xlabel(r"$f$ (Hz)")
plt.ylabel(r"$\mathscr{Y}$ (A/V)")
plt.grid(True)
legend.get_frame().set_facecolor('#00FFCC')

# Bottom subplot: Phase data comparison (if available)
plt.subplot(212)
plt.plot(xx, phi(xx, C0, R1, L1, C1), 'r-', label='phase model', linewidth=2)
if len(phase_data) > 0:
    plt.plot(phase_freq, phase_data, 'bo', label='phase data', markersize=4)
plt.xlabel(r"$f$ (Hz)")
plt.grid(True) 
plt.ylabel(r"$\phi$ (degrees)")
legend = plt.legend(loc='upper right', shadow=True, fontsize='large')
legend.get_frame().set_facecolor('#00FFCC')

# Save the plot
plt.tight_layout()
if plottype=='PNG' or plottype=='':
    plotname = infile.split('.')[0]+"trmodel"+loctime.replace(':','-')+'.PNG'
    plt.savefig(plotname, format='png', dpi=res)
else:
    plotname = infile.split('.')[0]+"trmodel"+loctime.replace(':','-')+'.eps'
    plt.savefig(plotname, format='eps', dpi=res)

print(f"Plot saved as: {plotname}")
plt.show()