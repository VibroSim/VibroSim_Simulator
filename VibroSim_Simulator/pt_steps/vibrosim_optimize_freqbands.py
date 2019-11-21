
import sys
import os
import numpy as np
from matplotlib import pyplot as pl

# Run this on output of modal analysis to interpret 
# the modal decay coefficients and plan a three-or four-segment
# frequency domain calculation that will be invertable 
# to a time-domain response. 


def run(dc_modalfreqs_href,
        dc_temporal_decay_divisor_float=50.0,
        dc_spectral_decay_divisor_float=10.0,):
    filename=dc_modalfreqs_href.getpath()


    min_freq = .5 # ignore normal modes below 0.5 Hz
    dt=1e-6 # 1 us timestep in final result

    
    (seg1_freqstart,seg1_freqstep,seg1_freqend,
     seg2_freqstart,seg2_freqstep,seg2_freqend,
     seg3_freqstart,seg3_freqstep,seg3_freqend,
     seg4_freqstart,seg4_freqstep,seg4_freqend
    ) = optimize_freqbands(filename,dc_temporal_decay_divisor_float,dc_spectral_decay_divisor_float,min_freq,dt,True)
    
    ret = {
        "dc:freqband_seg1_start": numericunitsv(seg1_freqstart,"Hz"),
        "dc:freqband_seg1_step": numericunitsv(seg1_freqstep,"Hz"),
        "dc:freqband_seg1_end": numericunitsv(seg1_freqend,"Hz"),
        "dc:freqband_seg2_start": numericunitsv(seg2_freqstart,"Hz"),
        "dc:freqband_seg2_step": numericunitsv(seg2_freqstep,"Hz"),
        "dc:freqband_seg2_end": numericunitsv(seg2_freqend,"Hz"),
        "dc:freqband_seg3_start": numericunitsv(seg3_freqstart,"Hz"),
        "dc:freqband_seg3_step": numericunitsv(seg3_freqstep,"Hz"),
        "dc:freqband_seg3_end": numericunitsv(seg3_freqend,"Hz"),
        }
    return ret
    
