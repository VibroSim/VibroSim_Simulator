
import sys
import os
import numpy as np
from matplotlib import pyplot as pl

# Run this on output of modal analysis to interpret 
# the modal decay coefficients and plan a three-or four-segment
# frequency domain calculation that will be invertable 
# to a time-domain response. 

from VibroSim_Simulator.read_spectrum import read_spectrum
from VibroSim_Simulator.optimize_freqbands import optimize_freqbands

def main(args=None):

    if args is None:
        args=sys.argv
        pass

    if len(args) < 2:
        print("Usage: %s <frequency table file>" % (args[0]))
        print("")
        print("Optimizes the frequency bands and step sizes for VibroSim simulation given modal analysis output.")
        print("")
        print("Accepts both COMSOL .txt output and .xls or .xlsx output with")
        print("\"Damped Frequency [Hz]\" and \"Logarithmic Decrement\" columns.")
        sys.exit(0)
        pass
    
    filename=args[1]

    temporal_decay_divisor = 50.0
    spectral_decay_divisor = 10.0

    min_freq = .5 # ignore normal modes below 0.5 Hz
    dt=1e-6 # 1 us timestep in final result


    complex_freqs_raw = read_spectrum(filename)
    (seg1_freqstart,seg1_freqstep,seg1_freqend,
     seg2_freqstart,seg2_freqstep,seg2_freqend,
     seg3_freqstart,seg3_freqstep,seg3_freqend,
     seg4_freqstart,seg4_freqstep,seg4_freqend
    ) = optimize_freqbands(complex_freqs_raw,temporal_decay_divisor,spectral_decay_divisor,min_freq,dt,True)
    
    pass
    
