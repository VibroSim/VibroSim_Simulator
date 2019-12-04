
import sys
import os
import numpy as np
from matplotlib import pyplot as pl

# Run this on output of modal analysis to interpret 
# the modal decay coefficients and plan a three-or four-segment
# frequency domain calculation that will be invertable 
# to a time-domain response. 

from VibroSim_Simulator.synthetic_spectrum import synthetic_spectrum
from VibroSim_Simulator.read_spectrum import read_spectrum

def main(args=None):

    if args is None:
        args=sys.argv
        pass

    if len(args) < 2:
        print("Usage: %s <frequency table file>" % (args[0]))
        print("")
        print("Plots a synthetic spectrum given modal analysis output.")
        print("")
        print("Accepts both COMSOL .txt output and .xls or .xlsx output with")
        print("\"Damped Frequency [Hz]\" and \"Logarithmic Decrement\" columns.")
        sys.exit(0)
        pass
        
    filename=args[1]

    complex_freqs_raw = read_spectrum(filename)
    synspec_fig = synthetic_spectrum(complex_freqs_raw)
    pl.show()
    
    pass
    
