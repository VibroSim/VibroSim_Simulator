import sys
import subprocess
import posixpath
import numpy as np
from matplotlib import pyplot as pl

from limatix.dc_value import hrefvalue as hrefv
from limatix.dc_value import numericunitsvalue as numericunitsv
from VibroSim_Simulator import read_comsol_probe_txt


def run(_xmldoc,_element, dc_modalfreqs_href):

    (metadata,fieldheaderdata,fieldrowdata) = read_comsol_probe_txt.read(dc_modalfreqs_href.getpath())
    complex_freqs_raw=fieldrowdata[0]['freq (1/s)'][1]



    df = 1.0 # frequency step in Hz
    fmax = 1.0e6

    # define frequency base
    frange = np.arange(fmax/df,dtype='d')*df

    synresp = np.zeros(frange.shape[0],dtype='D')
    
    for freqval in complex_freqs_raw:
        # Solution is np.exp(2*pi*i*freqval*t)u(t) or 
        # np.exp(2*pi*i*freqval.real*t - 2*pi*freqval.imag*t)*u(t)
        # Foruier transform pair: exp(-at)u(t) <==> 1/(a+2pi*i*f)
        # let a=-2pi*i*freqval=-2*pi*i*freqval.real+2*pi*freqval.imag
        # real(a) = 2*pi*freqval.imag
        # imag(a) = -2*pi*freqval.real
        # or a = real(a)+i*imag(a) = 2*pi*freqval.imag - 2*pi*freqval.real*i
        #      = (2*pi)*(freqval.imag -freqval.real*i)
        #      = -(2*pi)*i*freqval
        thisresp = 1.0/(-2.0*np.pi*(0+1j)*freqval + 2.0*np.pi*(0+1j)*frange)
        synresp += thisresp


        #pl.figure()
        #pl.clf()
        #pl.loglog(frange/1e3,np.abs(thisresp))
        #pl.xlabel('Frequency (kHz)')
        #pl.ylabel('Normalized response')
        #pl.grid()
        #pl.title('Synthetic spectrum for mode @ %f kHz' % (freqval.real/1e3))

        pass

    pl.figure()
    pl.clf()
    pl.loglog(frange,np.abs(synresp))
    maxabsresp = np.max(np.abs(synresp))
    pl.axis((10.0,100e3,maxabsresp/100.0,maxabsresp))
    pl.xlabel('Frequency (Hz)')
    pl.ylabel('Normalized response')
    pl.grid()
    pl.title('Synthetic spectrum from modal analysis')
    fig_href = hrefv(posixpath.split(dc_modalfreqs_href.get_bare_unquoted_filename())[0]+"_synthetic_spectrum.png",dc_modalfreqs_href.leafless())
    pl.savefig(fig_href.getpath(),dpi=300)
    
    ret = {
        "dc:modal_synthetic_spectrum": fig_href,
    }

    
    return ret
