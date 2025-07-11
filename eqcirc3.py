####################--
# File: eqcirc3.py
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


#from pylab import *
import sys
from scipy.optimize import leastsq
import numpy as np
from numpy import array
from numpy import sqrt
import matplotlib.pyplot as plt
import time #to allow time stamp on output
# Test for Python version:
cur_version = sys.version_info
# Initialize some lists:
ydat = []
x = []
xx = []
yy = []
zdat = []
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
    + 0.2e1*np.sqrt((fa ** 2 - f ** 2)**2/np.pi**4/fa**8*Ymin**4
    + 0.4e1*Ymin**2*Ymax**2/np.pi**4/fa**4))/0.4e1


def R1_i(Ymin, Ymax, fr, fa, C0):  # Motional resistance estimate
    return (-0.4e1*np.pi**2*fr**2*C0**2+Ymax**2)**(-0.1e1/0.2e1)


def L1_i(fr, fa, C0):  # Motional inductance estimate
    return 0.1e1 / np.pi ** 2 / (fa ** 2 - fr ** 2) / C0 / 0.4e1


def C1_i(fr, fa, C0):  # Motional capacitance estimate
    return (fa ** 2 / fr ** 2 - 1) * C0


def rez(z, ydat, f):  # Residual function
    return ydat - y(f, z)
#
# Comment out py3print() if you are using Python 2.x.x
#
def py3print():
    '''
        Author: D. S. Stutts
        4-22-2016
        This function uses the print() function
        according to the Python 3.x.x requirements.
        '''
    print( "Ymax = ", Ymax, " at fr = ", fr, "\n")
    print( "Ymin = ", Ymin, " at fa = ", fa, "\n")
    print( "fr = ", fr, "\n")
    print( "fa = ", fa, "\n")
    # Initial estimates:
    print("C0i = ", C0i,"\n")
    print("R1i = ", R1i,"\n")
    print("L1i = ", L1i,"\n")
    print("C1i = ", C1i,"\n")
    print("Qi = ", Qi,"\n")
    # Optimal estimates:
    print( "C0 = ", C0, "\n")
    print( "R1 = ", R1, "\n")
    print( "L1 = ", L1, "\n")
    print( "C1 = ", C1, "\n")
    print( "Q = ", Q, "\n")
    print( "k31 = ", k31, "\n")
    print("RMS Diviation = ", rmserr,"\n")

# Set the desired resolution:
res = 300# Use a larger value for PNG
#plottype = ''# Defaults to PNG
plottype = 'EPS'
# Input data file on command line:
infile = sys.argv[1]
data = open(infile, "r")  # get array out of input file
numline = 0
# Count the number of lines in the data file:
for line in data:
    numline +=1

# Calculate the number of magnitude data points:

data.seek(0) # Reset file pointer to the beginning
nummagpts = (numline - 1 - 26)/2
linecount = 0
# read the 21st through total lines from the data file
# and fill x,y lists with floating point numbers:
if cur_version[0]==3:# This is necesary due to the change in the type
    for line in data:# returned by the map function in Python 3.x.x.
        if linecount > 20 and linecount < nummagpts+21:# relative to 2.x.x.
            freqs = list(map(float, (line[0:31]).split()))
            impedances = list(map(float, (line[0:31]).split()))
            x.append(freqs[0])
            zdat.append(impedances[1])
        linecount += 1
else:
    for line in data:
        if linecount > 20 and linecount < nummagpts+21:
            x.append(map(float, (line[0:31]).split())[0])
            zdat.append(map(float, (line[0:31]).split())[1])
        linecount += 1

data.close()# close data file
xx = array(x)
zin = array(zdat)

yy = 1/zin

# Locate Ymax, Ymin, and initial guesses for fr and fa:
zmin = min(zdat)
zminidx = zdat.index(zmin)  # index of min impedance (max admittance)
zmax = max(zdat)
zmaxidx = zdat.index(zmax)  # index of max impedance (min admittance)
Ymax = max(yy)
Ymin = min(yy)
fr = x[zminidx]  # initial guess for resonance frequency
fa = x[zmaxidx]  # initial guess for antiresonance frequency
imax = len(x)

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

# Calculate RMS error: output is a tuple containing a dictionary
var = np.inner(output[2]['fvec'],output[2]['fvec'])
rmserr = sqrt(var/(len(output[2]['fvec'])-4))  # Standard error estimate
k31 = np.sqrt((fa**2-fr**2)/fa**2)# Electromechanical coupling coefficient
# Print the results to std out:
py3print()

# Calculate plot annotation positions:
delx = (fa-fr)/10.0
dely = (Ymax-Ymin)/20.0
noteymax = 0.65*Ymax
# Plot the model and the data:
plt.figure(figsize=(8,7),dpi=res)
plt.subplot(211)
plt.plot(xx, y(xx, coeffs), 'r-', label='model')
plt.plot(xx, yy, 'go', label='data')
#plt.errorbar(xx, yy, xerr=0, yerr=rmserr)# If accounting for x uncertainty
plt.errorbar(xx, y(xx, coeffs), yerr=rmserr, linestyle='None')
plt.annotate(r"$f_r$ = "+'{: 3.3e}'.format(fr),xy=(fa-delx,noteymax))
plt.annotate(r"$f_a$ = "+'{: 3.3e}'.format(fa),xy=(fa-delx,noteymax-1.4*dely))
plt.annotate(r"$C_0$ = "+'{: 3.3e}'.format(C0),xy=(fa-delx,noteymax-2.8*dely))
plt.annotate(r"$R_1$ = "+'{: 3.3e}'.format(R1),xy=(fa-delx,noteymax-4.4*dely))
plt.annotate(r"$L_1$ = "+'{: 3.3e}'.format(L1),xy=(fa-delx,noteymax-5.8*dely))
plt.annotate(r"$C_1$ = "+'{: 3.3e}'.format(C1),xy=(fa-delx,noteymax-7.2*dely))
plt.annotate('Q = '+'{: 3.3e}'.format(Q),xy=(fa-delx,noteymax-8.6*dely))
plt.annotate('RMS Dev. = '+'{: 3.2e}'.format(rmserr),xy=(fa-delx,noteymax-10*dely))
plt.annotate('$k_{31}$ = '+'{: 3.2e}'.format(k31),xy=(fa-delx,noteymax-11.4*dely))

# Add date and time in plot title:
loctime = time.asctime(time.localtime(time.time()))
plt.suptitle('Data File ='+infile+':  '+loctime)
#plt.annotate('Date and time = '+localtime,xy=(fa,noteymax+dely))
print ("Date and Time =", loctime, "\n")


legend = plt.legend(loc='upper right', shadow=True, fontsize='large')
plt.xlabel(r"$f$ (Hz)")
plt.ylabel(r"$\mathscr{Y}$ (A/V)")
plt.grid(True)
# Put a nice background color on the legend:
legend.get_frame().set_facecolor('#00FFCC')
plt.subplot(212)
plt.plot(xx, phi(xx, C0, R1, L1, C1), 'r-', label='model')
plt.xlabel(r"$f$ (Hz)")
plt.grid(True)
plt.ylabel(r"$\phi$ (degrees)")
legend = plt.legend(loc='upper right', shadow=True, fontsize='large')
legend.get_frame().set_facecolor('#00FFCC')
if plottype=='PNG' or plottype=='':# Default to PNG
    # Save plot as PNG:
    plotname = infile.split('.')[0]+"trmodel"+loctime.replace(':','-')+'.PNG'
    plt.savefig(plotname,format='png', dpi=res)
else:# Save plot as EPS:
    plotname = infile.split('.')[0]+"trmodel"+loctime.replace(':','-') + '.eps'
    plt.savefig(plotname,format='eps', dpi=res)

plt.show()

