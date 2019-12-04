import sys
import subprocess
import posixpath
import numpy as np
from matplotlib import pyplot as pl

from limatix.dc_value import hrefvalue as hrefv
from limatix.dc_value import numericunitsvalue as numericunitsv

from VibroSim_Simulator.synthetic_spectrum import synthetic_spectrum
from VibroSim_Simulator.read_spectrum import read_spectrum

def run(_xmldoc,_element, dc_modalfreqs_href):

    complex_freqs_raw = read_spectrum(dc_modalfreqs_href.getpath())

    synspec_fig = synthetic_spectrum(complex_freqs_raw)

    # Save figure to disk
    pl.figure(synspec_fig.fignum)
    fig_href = hrefv(posixpath.splitext(dc_modalfreqs_href.get_bare_unquoted_filename())[0]+"_synthetic_spectrum.png",dc_modalfreqs_href.leafless())
    pl.savefig(fig_href.getpath(),dpi=300)
    
    ret = {
        "dc:modal_synthetic_spectrum": fig_href,
    }

    
    return ret
