import sys
import subprocess
import numpy as np

from limatix.dc_value import numericunitsvalue as numericunitsv
from VibroSim_Simulator import format_modes
from VibroSim_Simulator import enter_frequency


def run(_xmldoc,_element, dc_modalfreqs_href):

    (metadata,fieldheaderdata,fieldrowdata) = read_comsol_probe_txt.read(dc_modalfreqs_href)
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
        synresp += 1.0/(-2.0*np.pi*(0+1j)*freqval + 2.0*np.pi*(0+1j)+frange)
        pass

    pl.figure()
    pl.clf()
    pl.loglog(frange/1e3,synresp)
    pl.xlabel('Frequency (kHz)')
    pl.ylabel('Normalized response')
    pl.grid()
    pl.title('Synthetic spectrum from modal analysis')
    fig_href = hrefv(posixpath.split(dc_modalfreqs_href.get_bare_unquoted_filename)[0]+"_synthetic_spectrum.png",dc_modalfreqs_href.leafless())
    pl.savefig(fig_href.getpath(),dpi=300)
    
    
    return {
        "dc:modal_synthetic_spectrum": fig_href,
    }
