import sys
import os
import os.path

try:
    # py2.x
    from urllib import pathname2url
    from urllib import url2pathname
    from urllib import quote
    from urllib import unquote
    pass
except ImportError:
    # py3.x
    from urllib.request import pathname2url
    from urllib.request import url2pathname
    from urllib.parse import quote
    from urllib.parse import unquote
    pass

import numpy as np
import pandas as pd

from limatix.dc_value import hrefvalue as hrefv

from VibroSim_WelderModel.contact_model import write_motiontable


# vibrosim_calc_singlefrequency_motion.py:
# Processtrak for evaluating motion corresponding
# to single-frequency excitation. 
#
# This motion can then be used by vibrosim_calc_heating_welder.py
# which should give approximately the same heating as 
# vibrosim_calc_heating_singlefrequency.py. It makes 
# a good cross-check

def run(_xmldoc,_element, 
        dc_dest_href,
        dc_measident_str,
        dc_exc_t0_numericunits,
        dc_exc_t1_numericunits,
        dc_exc_t2_numericunits,
        dc_exc_t3_numericunits,
        dc_exc_t4_numericunits,
        dc_excitation_frequency_numericunits,
        dc_harmonicburst_normalstress_complex, 
        dc_harmonicburst_shearstress_complex):

    t0 = dc_exc_t0_numericunits.value("s")
    t1 = dc_exc_t1_numericunits.value("s")
    t2 = dc_exc_t2_numericunits.value("s")
    t3 = dc_exc_t3_numericunits.value("s")
    t4 = dc_exc_t4_numericunits.value("s")

    freq = dc_excitation_frequency_numericunits.value("Hz")

    dt=1e-6

    trange = np.arange(dc_exc_t4_numericunits.value("s")/dt)*dt
    
    heatingenvelope = np.zeros(trange.shape[0],dtype='d')
    heatingenvelope[ (trange >= t0) & (trange < t1) ] = (1.0-np.cos( np.pi*(trange[ (trange >= t0) & (trange < t1) ]-t0)/(t1-t0)))/2.0
    heatingenvelope[ (trange >= t1) & (trange < t2) ] = 1.0
    heatingenvelope[ (trange >= t2) & (trange < t3) ] = (1.0+np.cos( np.pi*(trange[ (trange >= t2) & (trange < t3) ]-t2)/(t3-t2)))/2.0

    crackcenternormalstress = np.real(heatingenvelope*dc_harmonicburst_normalstress_complex*np.exp((0+1j)*2.0*np.pi*freq*(trange-t0)))
    crackcentershearstress = np.real(heatingenvelope*dc_harmonicburst_shearstress_complex*np.exp((0+1j)*2.0*np.pi*freq*(trange-t0)))

    motiontable = pd.DataFrame(index=pd.Float64Index(data=np.arange(trange.shape[0],dtype='d')*dt,dtype='d',name="Time(s)"))
    motiontable.insert(len(motiontable.columns),"specimen_crackcenternormalstress(Pa)",crackcenternormalstress)
    motiontable.insert(len(motiontable.columns),"specimen_crackcentershearstress(Pa)",crackcentershearstress)


    # Save motiontable CSV and add to return dictionary
    motiontable_href = hrefv(quote("%s_motiontable.csv.bz2" % (dc_measident_str)),dc_dest_href)
    write_motiontable(motiontable,motiontable_href.getpath())

    ret = {
        "dc:motion": motiontable_href
    }

    return ret

