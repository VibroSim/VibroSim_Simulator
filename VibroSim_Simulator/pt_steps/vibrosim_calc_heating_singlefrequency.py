import sys
import numpy as np
import sys
import os
import os.path
import tempfile
import scipy 
import numpy as np
import pandas as pd
from matplotlib import pyplot as pl

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

from limatix.dc_value import hrefvalue as hrefv
from limatix.dc_value import numericunitsvalue as numericunitsv


from crackclosuresim2 import inverse_closure
from crackclosuresim2 import crackopening_from_tensile_closure
from crackclosuresim2 import solve_normalstress
from crackclosuresim2 import load_closurestress


from angled_friction_model.angled_friction_model import angled_friction_model
from angled_friction_model.angled_friction_model import integrate_power

from limatix.dc_value import numericunitsvalue as numericunitsv

from VibroSim_Simulator.calc_heating import calc_heating_singlefrequency


def run(_xmldoc,_element, 
        dc_dest_href,
        dc_measident_str,
        dc_friction_coefficient_float,
        dc_msqrtR_numericunits,
        dc_staticload_numericunits,
        dc_exclength_numericunits,
        dc_tortuosity_numericunits,
        dc_numdraws_int,
        dc_YoungsModulus_numericunits,
        dc_PoissonsRatio_float,
        dc_YieldStrength_numericunits,
        dc_Density_numericunits,
        dc_closurestate_side1_href,
        dc_closurestate_side2_href,
        dc_a_side1_numericunits,
        dc_a_side2_numericunits,
        dc_crack_model_normal_str,
        dc_crack_model_shear_str,
        dc_crack_model_shear_factor_float, # shear sensitivity factor (nominally 1.0)
        dc_exc_t0_numericunits,
        dc_exc_t1_numericunits,
        dc_exc_t2_numericunits,
        dc_exc_t3_numericunits,
        dc_exc_t4_numericunits,
        dc_excitation_frequency_numericunits,
        dc_harmonicburst_normalstress_complex, 
        dc_harmonicburst_shearstress_complex):

    verbose=False
    doplots=False


    ## Manually extract dc_heatingdata_href so we can do it with no provenance
    #dc_heatingdata_el = _xmldoc.xpathsinglecontext(_element,"dc:heatingdata")
    #dc_heatingdata_href=hrefv.fromxml(_xmldoc,dc_heatingdata_el)

    dc_heatingdata_href = hrefv(quote(dc_measident_str+"_heatingdata.txt"),dc_dest_href)


    (xrange_side1,
     x_bnd_side1,
     dx_side1,
     a_side1_verify,
     closure_stress_side1,
     crack_opening_side1) = load_closurestress(dc_closurestate_side1_href.getpath())

    if a_side1_verify is not None: # crack lengths should match
        assert((a_side1_verify-dc_a_side1_numericunits.value("m"))/dc_a_side1_numericunits.value("m") < 1e-2)
        pass

    (xrange_side2,
     x_bnd_side2,
     dx_side2,
     a_side2_verify,
     closure_stress_side2,
     crack_opening_side2) = load_closurestress(dc_closurestate_side2_href.getpath())

    if a_side2_verify is not None: # crack lengths should match
        assert((a_side2_verify-dc_a_side2_numericunits.value("m"))/dc_a_side2_numericunits.value("m") < 1e-2)
        pass

    #closurestate_side1_dataframe = pd.read_csv(dc_closurestate_side1_href.getpath(),index_col=0)
    #closurestate_side2_dataframe = pd.read_csv(dc_closurestate_side2_href.getpath(),index_col=0)


    #xrange_side1 = np.array(closurestate_side1_dataframe.index)
    
    #closure_stress_side1 = np.array(closurestate_side1_dataframe["Closure stress (Pa)"])
    #crack_opening_side1 = np.array(closurestate_side1_dataframe["Crack opening (m)"])


    #xrange_side2 = np.array(closurestate_side2_dataframe.index)
    #    
    #closure_stress_side2 = np.array(closurestate_side2_dataframe["Closure stress (Pa)"])
    #crack_opening_side2 = np.array(closurestate_side2_dataframe["Crack opening (m)"])

    shorter_xrange_len = min(xrange_side1.shape[0],xrange_side2.shape[0])

    if np.any(xrange_side1[:shorter_xrange_len] != xrange_side2[:shorter_xrange_len]):
        raise ValueError("Crack radius positions from %s and %s do not match!" % (dc_closurestate_side1_href.humanurl(),dc_closurestate_side2_href.humanurl()))

    
    # xrange should be the longer of the two possibilities
    if xrange_side1.shape[0]==shorter_xrange_len:
        xrange = xrange_side2
        x_bnd = x_bnd_side2
        pass
    else:
        xrange=xrange_side1
        x_bnd = x_bnd_side1
        pass
    
    xstep = xrange[1]-xrange[0]

    # Extend fields out appropriately for our xrange
    closure_stress_side1_ext = np.concatenate((closure_stress_side1,np.ones(xrange.shape[0]-closure_stress_side1.shape[0],dtype='d')*closure_stress_side1[-1]))
    closure_stress_side2_ext = np.concatenate((closure_stress_side2,np.ones(xrange.shape[0]-closure_stress_side2.shape[0],dtype='d')*closure_stress_side2[-1]))

    crack_opening_side1_ext = np.concatenate((crack_opening_side1,np.zeros(xrange.shape[0]-closure_stress_side1.shape[0],dtype='d')))
    crack_opening_side2_ext = np.concatenate((crack_opening_side2,np.zeros(xrange.shape[0]-closure_stress_side2.shape[0],dtype='d')))


    # can scriptify(calc_heating_singlefrequency) to help with debugging
    (heatpower_fig,
     power_per_m2_side1,
     power_per_m2_side2,
     totalpower) = calc_heating_singlefrequency(dc_friction_coefficient_float,
                                                dc_msqrtR_numericunits.value("m^-1.5"),
                                                dc_staticload_numericunits.value("Pa"),
                                                dc_tortuosity_numericunits.value("degrees")*np.pi/180.0, # angular_stddev
                                                dc_numdraws_int,
                                                
                                                dc_YoungsModulus_numericunits.value("Pa"),
                                                dc_PoissonsRatio_float,
                                                dc_YieldStrength_numericunits.value("Pa"),
                                                dc_Density_numericunits.value("kg/m^3"),
                                                x_bnd,xrange,xstep, # base of crack radius values
                                                closure_stress_side1_ext, # side1 closure state
                                                crack_opening_side1_ext,
                                                dc_a_side1_numericunits.value("m"),                        
                                                closure_stress_side2_ext, # side2 clsoure state
                                                crack_opening_side2_ext,
                                                dc_a_side2_numericunits.value("m"),                        
                                                dc_crack_model_normal_str,
                                                dc_crack_model_shear_str,
                                                dc_crack_model_shear_factor_float, # shear sensitivity factor (nominally 1.0)
                                                dc_excitation_frequency_numericunits.value("Hz"),
                                                                                                
                                                dc_exc_t0_numericunits.value('s'),
                                                dc_exc_t1_numericunits.value('s'),
                                                dc_exc_t2_numericunits.value('s'),
                                                dc_exc_t3_numericunits.value('s'),
                                                dc_exc_t4_numericunits.value('s'),
                                                dc_harmonicburst_normalstress_complex,
                                                dc_harmonicburst_shearstress_complex,
                                                dc_heatingdata_href.getpath())
    
    
    
    pl.figure(heatpower_fig.number)
    heatpower_href = hrefv(quote(dc_measident_str+"_heatpower.png"),dc_dest_href)
    pl.savefig(heatpower_href.getpath(),dpi=300)
     
    

    ret = { 
        "dc:heatpower": heatpower_href,
        "dc:heatingdata": dc_heatingdata_href,
        "dc:heatingtotalpower": numericunitsv(totalpower,"W"),
    }

    return ret
