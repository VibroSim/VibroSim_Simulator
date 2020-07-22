import sys
import os
import numpy as np

from . import read_comsol_probe_txt

def print_comsol_modes(filename):
    (metadata,fieldheaderdata,fieldrowdata) = read_comsol_probe_txt.read(filename)

    if "freq (1/s)" in fieldrowdata[0]:
        complex_freqs_raw=fieldrowdata[0]['freq (1/s)'][1]
        pass
    elif "freq (Hz)" in fieldrowdata[0]: # COMSOL 5.5 gives freq in Hz instead of 1/s as column label
        complex_freqs_raw=fieldrowdata[0]['freq (Hz)'][1]
        pass
    else: 
        raise ValueError("Did not find field row data for frequencies (candidates=%s)" % (str(list(fieldrowdata[0].keys()))))

    print("Frequency (Hz)   1/e Decay time (s)")
    print("-----------------------------------")
    print("\n".join(["%15.2f  %15.2f" % (freqval.real,1.0/(2.0*np.pi*freqval.imag)) for freqval in complex_freqs_raw])) 
    print("-----------------------------------")
    pass
    
