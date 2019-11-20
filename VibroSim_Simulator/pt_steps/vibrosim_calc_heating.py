import sys
import numpy as np
import sys
import os
import os.path
import tempfile
import scipy 
import numpy as np
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

from crackclosuresim2 import ModeI_throughcrack_CODformula
from crackclosuresim2 import Tada_ModeI_CircularCrack_along_midline
from crackclosuresim2 import ModeII_throughcrack_CSDformula
from crackclosuresim2.fabrikant import Fabrikant_ModeII_CircularCrack_along_midline

from angled_friction_model.angled_friction_model import angled_friction_model
from angled_friction_model.angled_friction_model import integrate_power



from limatix.dc_value import numericunitsvalue as numericunitsv


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
	dc_reff_side1_array, # NOTE: arrayvalue class is not unit-aware!
	dc_seff_side1_array, 
	dc_reff_side2_array,
	dc_seff_side2_array,
	dc_exc_t0_numericunits,
	dc_exc_t1_numericunits,
	dc_exc_t2_numericunits,
	dc_exc_t3_numericunits,
	dc_excitation_frequency_numericunits,
        dc_harmonicburst_normalstrain_complex, 
        dc_harmonicburst_shearstrain_complex):

    verbose=False
    doplots=False


    # Manually extract dc_heatingdata_href so we can do it with no provenance
    dc_heatingdata_el = _xmldoc.xpathsinglecontext(_element,"dc:heatingdata")
    dc_heatingdata_href=hrefv.fromxml(_xmldoc,dc_heatingdata_el)

    t=np.linspace(dc_exc_t0_numericunits.value('s'),dc_exc_t3_numericunits.value('s'),100)
    
    friction_coefficient = dc_friction_coefficient_float
    msqrtR = dc_msqrtR_numericunits.value("m^-1.5")
    static_load = dc_staticload_numericunits.value("Pa")

    plane_stress_modulus = dc_YoungsModulus_numericunits.value("Pa")/(1.0-dc_PoissonsRatio_float**2.0)
    shear_modulus = dc_YoungsModulus_numericunits.value("Pa")/(2.0*(1.0+dc_PoissonsRatio_float))

    # Note: We are ignoring the phases of the normal and shear strains
    vib_normal_stress_ampl = plane_stress_modulus * abs(dc_harmonicburst_normalstrain_complex)
    vib_shear_stress_ampl = shear_modulus * abs(dc_harmonicburst_shearstrain_complex)

    
    angular_stddev = dc_tortuosity_numericunits.value("degrees")*np.pi/180.0

    sigma_yield = dc_YieldStrength_numericunits.value("Pascals")
    tau_yield = sigma_yield/2.0
    
    #crack_model_normal = ModeI throughcrack_CODformula(E)
    crack_model_normal = Tada_ModeI_CircularCrack_along_midline(dc_YoungsModulus_numericunits.value("Pa"),dc_PoissonsRatio_float)
    #crack_model_shear = ModeII_throughcrack_CSDformula(E,nu)
    crack_model_shear = Fabrikant_ModeII_CircularCrack_along_midline(dc_YoungsModulus_numericunits.value("Pa"),dc_PoissonsRatio_float)


    reff_side1 = dc_reff_side1_array.value()
    seff_side1 = dc_seff_side1_array.value()

    reff_side2 = dc_reff_side2_array.value()
    seff_side2 = dc_seff_side2_array.value()


    # Fully open crack lengths for left and right side
    aside1 = np.max(reff_side1) 
    aside2 = np.max(reff_side2)

    # Desired approximate step size for calculations
    approximate_xstep=25e-6 # 25um

    num_boundary_steps=np.floor((max(aside1,aside2)+approximate_xstep)/approximate_xstep)
    xmax = num_boundary_steps*approximate_xstep  # Maximum position from center to calculate to;
    # should exceed half-crack lengths 

    numsteps = num_boundary_steps-1
    xstep = (xmax)/(numsteps) # Actual step size so that xmax is a perfect multiple of this number

    x_bnd = xstep*np.arange(num_boundary_steps) # Position of element boundaries
    xrange = (x_bnd[1:] + x_bnd[:-1])/2.0 # Position of element centers

    # Determine closure stress field from observed crack length data
    closure_stress_side1=inverse_closure(reff_side1,seff_side1,xrange,x_bnd,xstep,aside1,sigma_yield,crack_model_normal,verbose=verbose)
    
    closure_stress_side2=inverse_closure(reff_side2,seff_side2,xrange,x_bnd,xstep,aside2,sigma_yield,crack_model_normal,verbose=verbose)
    
    
    # Evaluate initial crack opening gaps from extrapolated tensile closure field
    crack_initial_opening_side1 = crackopening_from_tensile_closure(xrange,x_bnd,closure_stress_side1,xstep,aside1,sigma_yield,crack_model_normal)
    
    crack_initial_opening_side2 = crackopening_from_tensile_closure(xrange,x_bnd,closure_stress_side2,xstep,aside2,sigma_yield,crack_model_normal)
    
    # Plot the evaluated closure state (side1)
    pl.figure()
    pl.plot(xrange[xrange < aside1]*1e3,closure_stress_side1[xrange < aside1]/1e6,'-',
            reff_side1*1e3,seff_side1/1e6,'x')
    for observcnt in range(len(reff_side1)):        
        (effective_length, sigma, tensile_displ, dsigmaext_dxt) = solve_normalstress(xrange,x_bnd,closure_stress_side1,xstep,seff_side1[observcnt],aside1,sigma_yield,crack_model_normal)
	pl.plot(effective_length*1e3,seff_side1[observcnt]/1e6,'.')
	pass
    pl.grid(True)
    pl.legend(('Closure stress field','Observed crack tip posn','Recon. crack tip posn'),loc="best")
    pl.xlabel('Radius from crack center (mm)')
    pl.ylabel('Stress (MPa)')
    pl.title('Crack closure state (side1)')
    
    closure_state_side1_href = hrefv(quote(dc_measident_str+"_closurestate_side1.png"),dc_dest_href)
    pl.savefig(closure_state_side1_href.getpath(),dpi=300)

    # Plot the evaluated closure state (side2)
    pl.figure()
    pl.plot(xrange[xrange < aside2]*1e3,closure_stress_side2[xrange < aside2]/1e6,'-',
            reff_side2*1e3,seff_side2/1e6,'x')
    for observcnt in range(len(reff_side2)):        
        (effective_length, sigma, tensile_displ, dsigmaext_dxt) = solve_normalstress(xrange,x_bnd,closure_stress_side2,xstep,seff_side2[observcnt],aside2,sigma_yield,crack_model_normal)
	pl.plot(effective_length*1e3,seff_side2[observcnt]/1e6,'.')
	pass
    pl.grid(True)
    pl.legend(('Closure stress field','Observed crack tip posn','Recon. crack tip posn'),loc="best")
    pl.xlabel('Radius from crack center (mm)')
    pl.ylabel('Stress (MPa)')
    pl.title('Crack closure state (side 2)')
    
    closure_state_side2_href = hrefv(quote(dc_measident_str+"_closurestate_side2.png"),dc_dest_href)
    pl.savefig(closure_state_side2_href.getpath(),dpi=300)



    # Now calculate crack heating under the given conditions
    (power_per_m2_side1,
     power_per_m2_stddev_side1,
     vibration_ampl_side1) = angled_friction_model(x_bnd,xrange,xstep,
                                                   dc_numdraws_int,
						   dc_YoungsModulus_numericunits.value("Pa"),
						   dc_PoissonsRatio_float,
						   sigma_yield,tau_yield,
						   friction_coefficient,
						   closure_stress_side1,
						   crack_initial_opening_side1,
						   angular_stddev,
						   aside1,
						   static_load,
						   vib_normal_stress_ampl,
						   vib_shear_stress_ampl,
						   dc_excitation_frequency_numericunits.value("Hz"),
						   crack_model_normal,
						   crack_model_shear,
						   1.0,
						   msqrtR,
						   verbose,
						   doplots)
						  
    (power_per_m2_side2,
     power_per_m2_stddev_side2,
     vibration_ampl_side2) = angled_friction_model(x_bnd,xrange,xstep,
                                                   dc_numdraws_int,
						   dc_YoungsModulus_numericunits.value("Pa"),
						   dc_PoissonsRatio_float,
						   sigma_yield,tau_yield,
						   friction_coefficient,
						   closure_stress_side2,
						   crack_initial_opening_side2,
						   angular_stddev,
						   aside2,
						   static_load,
						   vib_normal_stress_ampl,
						   vib_shear_stress_ampl,
						   dc_excitation_frequency_numericunits.value("Hz"),
						   crack_model_normal,
						   crack_model_shear,
						   1.0,
						   msqrtR,
						   verbose,
						   doplots)

    (totalpower_side1, totalpower_stddev_side1) = integrate_power(xrange,power_per_m2_side1,power_per_m2_stddev_side1)
    (totalpower_side2, totalpower_stddev_side2) = integrate_power(xrange,power_per_m2_side2,power_per_m2_stddev_side2)

    totalpower=totalpower_side1 + totalpower_side2
     
    pl.figure()
    pl.clf()
    pl.plot(-xrange*1e3,power_per_m2_side1/1.e3,'-',
            xrange*1e3,power_per_m2_side2/1.e3,'-',)
    pl.grid()
    pl.xlabel('Radius from center (mm)')
    pl.ylabel('Heating power (kW/m^2)')
    pl.title('Crack power deposition')

    heatpower_href = hrefv(quote(dc_measident_str+"_heatpower.png"),dc_dest_href)
    pl.savefig(heatpower_href.getpath(),dpi=300)

    heatingfh = open(dc_heatingdata_href.getpath(),"w")
    heatingfh.write("% t(s) \t r(m) \t side1_heating(W/m^2) \t side2_heating(W/m^2)\n")

    for tidx in range(t.shape[0]):
        for xidx in range(xrange.shape[0]):
            heatingfh.write("%.8e\t%.8e\t%.8e\t%.8e\n" % (t[tidx],xrange[xidx],power_per_m2_side1[xidx],power_per_m2_side2[xidx]))
            pass
        pass
    heatingfh.close()
            
    

    return { 
        "dc:closureplot_side1": closure_state_side1_href,
        "dc:closureplot_side2": closure_state_side2_href,
        "dc:heatpower": heatpower_href,
        "dc:heatingdata": dc_heatingdata_href,
        "dc:heatingtotalpower": numericunitsv(totalpower,"W"),
    }
