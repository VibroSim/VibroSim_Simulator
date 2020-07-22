import sys
import os
import os.path
import numpy as np
from matplotlib import pyplot as pl


def read_spectrum(filename):
    if filename.lower().endswith(".txt"):
        # .txt extension... assume COMSOL
        from VibroSim_Simulator import read_comsol_probe_txt
        (metadata,fieldheaderdata,fieldrowdata) = read_comsol_probe_txt.read(filename)
        if "freq (1/s)" in fieldrowdata[0]:
            complex_freqs_raw=fieldrowdata[0]['freq (1/s)'][1]
            pass
        elif "freq (Hz)" in fieldrowdata[0]: # COMSOL 5.5 gives freq in Hz instead of 1/s as column label
            complex_freqs_raw=fieldrowdata[0]['freq (Hz)'][1]
            pass
        else: 
            raise ValueError("Did not find field row data for frequencies (candidates=%s)" % (str(list(fieldrowdata[0].keys()))))
        pass
    elif filename.lower().endswith(".xls") or filename.lower().endswith(".xlsx"):
        # Excel... assume ANSYS
        import pandas as pd

        
        dataframe=pd.read_excel(filename)

        # Construct complex frequency from damped frequency and
        # logarithmic decrement
        freqs_realpart = np.array(dataframe["Damped Frequency [Hz]"])
        freqs_logdec = np.array(dataframe["Logarithmic Decrement"])
        # logarithmic decrement is ln(x2/x1) defined as the
        # logarithm of the ratio amplitudes of successive peaks
        # of the wave
        
        # So if we have the wave exp(i*2*pi*(f_real + i*f_imag)*t)
        # over one cycle starting at t=0 the 2*pi*f_real*t
        # term in the exponent increases by exactly 2 pi.
        # equivalently f_real*t = 1 or t=1/f_real.
        # Plugging this in, the wave value after one cycle becomes
        # exp(i*2*pi*(1 + i*f_imag/f_real))
        # with the decrement coming from the real part
        # exp(-2*pi*f_imag/f_real)
        #
        # ... which by definition the parameter of the exponent
        # must be the logarithmic decrement.
        #
        # so we have -2*pi*f_imag/f_real = logdec
        # or f_imag = -logdec*f_real/(2*pi)
        freqs_imagpart = -freqs_logdec*freqs_realpart/(2.0*np.pi)
        complex_freqs_raw = freqs_realpart + (0+1j)*freqs_imagpart
        pass
    else:
        raise ValueError("Unknown filename extension %s" % (os.path.splitext(filename)[1]))

    return complex_freqs_raw
