import sys
import os
import numpy as np

from . import read_comsol_probe_txt

def print_comsol_modes(filename):
    (metadata,fieldheaderdata,fieldrowdata) = read_comsol_probe_txt.read(filename)
    complex_freqs_raw=fieldrowdata[0]['freq (1/s)'][1]
    print("Frequency (Hz)   1/e Decay time (s)")
    print("-----------------------------------")
    print("\n".join(["%15.2f  %15.2f" % (freqval.real,1.0/(2.0*np.pi*freqval.imag)) for freqval in complex_freqs_raw])) 
    print("-----------------------------------")
    pass
    
