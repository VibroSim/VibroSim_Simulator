import sys
import numpy as np
import sys
import os
import os.path
import tempfile
import inspect
import scipy 
import scipy.interpolate
from scipy.interpolate import splrep,splev
import scipy.signal
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

from crackclosuresim2 import ModeI_throughcrack_CODformula
from crackclosuresim2 import Tada_ModeI_CircularCrack_along_midline
from crackclosuresim2 import ModeII_throughcrack_CSDformula
from crackclosuresim2.fabrikant import Fabrikant_ModeII_CircularCrack_along_midline

from angled_friction_model.angled_friction_model import angled_friction_model
from angled_friction_model.angled_friction_model import integrate_power



from limatix.dc_value import numericunitsvalue as numericunitsv


def fit_surrogate(strainvals,heatingvals):
    # Fit surrogate with splrep()
    # Allow two interior knots and a cubic spline fit
    # Fit over axis 0 of heatingvals (corresponds to axis of strainvals); iterate over axis 1 of heatingvals

    min_sv = np.min(strainvals)
    max_sv = np.max(strainvals)

    interior_knots=(
        min_sv + (max_sv-min_sv)/3.0,
        min_sv + (max_sv-min_sv)*2.0/3.0,
    )
    

    tck = np.zeros(heatingvals.shape[1],dtype='O')
    
    for xidx in range(heatingvals.shape[1]):
        tck[xidx] = scipy.interpolate.splrep(strainvals,heatingvals[:,xidx],
                                             k=3,
                                             task=-1,
                                             t=interior_knots)
        pass
    return tck
        

def filterfine(finedata,finetrange,coarsetrange):
    # Filter and coarsen along axis 0 of finedata
    
    fine_dt = finetrange[1]-finetrange[0]
    coarse_dt = coarsetrange[1]-coarsetrange[0]

    coarsen_factor=int(coarse_dt/fine_dt)

    # running average filter
    filter_npoints = coarsen_factor*2-1 # odd

    # Only recent scipy's have axes parameter to vectorize fftconvolve. 
    # otherwise we have to iterate
    if "axes" in inspect.getargspec(scipy.signal.fftconvolve).args:
        filtered_fine_center = scipy.signal.fftconvolve(finedata,np.ones(filter_npoints,dtype='d')/filter_npoints,mode='valid',axes=0) # length m-n+1
        pass
    else:
        filtered_fine_center=np.zeros((finedata.shape[0]-filter_npoints+1,finedata.shape[1]),dtype='d')
        for ax1pos in range(finedata.shape[1]):
            filtered_fine_center[:,ax1pos] = scipy.signal.fftconvolve(finedata[:,ax1pos],np.ones(filter_npoints,dtype='d')/filter_npoints,mode='valid') # length m-n+1
            pass
        pass

    filtered_fine_leftedge = np.zeros(((filter_npoints-1)//2,finedata.shape[1]),dtype='d')  # length (n-1)//2
    filtered_fine_rightedge = np.zeros(((filter_npoints-1)//2,finedata.shape[1]),dtype='d')  # length (n-1)//2

    for edgeidx in range((filter_npoints-1)//2):
        filtered_fine_leftedge[edgeidx:] = np.mean(finedata[:(edgeidx*2+1),:],axis=0)
        filtered_fine_rightedge[(-edgeidx-1),:] = np.mean(finedata[(-edgeidx*2-1):,:],axis=0)
        pass
    filtered_fine = np.concatenate((filtered_fine_leftedge,filtered_fine_center,filtered_fine_rightedge),axis=0)

    filtered_coarse=scipy.interpolate.interp1d(finetrange,filtered_fine,fill_value=0.0,bounds_error=False,assume_sorted=True,axis=0)(coarsetrange)
    
    # total length: m-n+1+(n-1)/2+(n-1)/2 = m-n+1+n-1 = m 


    return filtered_coarse

    
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
        dc_exc_t4_numericunits,
        dc_motion_href,
):

    verbose=False
    doplots=False


    # Manually extract dc_heatingdata_href so we can do it with no provenance
    dc_heatingdata_el = _xmldoc.xpathsinglecontext(_element,"dc:heatingdata")
    dc_heatingdata_href=hrefv.fromxml(_xmldoc,dc_heatingdata_el)

    t=np.linspace(dc_exc_t0_numericunits.value('s'),dc_exc_t4_numericunits.value('s'),100)
    
    friction_coefficient = dc_friction_coefficient_float
    msqrtR = dc_msqrtR_numericunits.value("m^-1.5")
    static_load = dc_staticload_numericunits.value("Pa")

    plane_stress_modulus = dc_YoungsModulus_numericunits.value("Pa")/(1.0-dc_PoissonsRatio_float**2.0)
    shear_modulus = dc_YoungsModulus_numericunits.value("Pa")/(2.0*(1.0+dc_PoissonsRatio_float))

    # Load motion file
    motiontable = pd.read_csv(dc_motion_href.getpath(),index_col=0)

    
    trange = np.array(motiontable.index) # ["Time(s)"]
    dt=trange[1]-trange[0]
    

    # ***!!! Should crosscheck with single frequency result by evaluating a motion table from single frequency data
    # and running that through this !!!***
    crack_normalstrain = np.array(motiontable["specimen_crackcenternormalstrain"])
    crack_shearstrain = np.array(motiontable["specimen_crackcentershearstrain"])
    
    ## pad to power-of-two size
    #n = 2.0**np.ceil(np.log(trange_orig.shape[0]*1.0)/np.log(2.0)) # next power of 2
    #
    #trange=trange_orig[0]+np.arange(n,dtype='d')*dt
    # 
    #crack_normalstrain=np.concatenate((crack_normalstrain_orig,np.zeros(n-trange_orig.shape[0],dtype='d')))
    #crack_shearstrain=np.concatenate((crack_shearstrain_orig,np.zeros(n-trange_orig.shape[0],dtype='d')))
    #
    #df=1.0/(n*dt)
    #frange=np.arange(n,dtype='d')*df
    #
    #crack_normalstrain_fft = np.fft.fft(crack_normalstrain)
    #crack_normalstrain_fft = np.fft.fft(crack_shearstrain)

    # To simplify matters we treat the normal and shear vibrations as uncoupled and separate. 


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



    # Find extrema
    normalstrain_relative_maxima = (crack_normalstrain[1:-1] > crack_normalstrain[2:]) & (crack_normalstrain[1:-1] > crack_normalstrain[:-2])
    normalstrain_relative_minima = (crack_normalstrain[1:-1] < crack_normalstrain[2:]) & (crack_normalstrain[1:-1] < crack_normalstrain[:-2])
    normalstrain_extrema = normalstrain_relative_minima | normalstrain_relative_maxima
    normalstrain_extrema_indices = np.concatenate(((0,),np.where(normalstrain_extrema)[0]+1,(crack_normalstrain.shape[0]-1,)))
    
    normalstrain_num_extrema=normalstrain_extrema_indices.shape[0]-1

    normalstrain_segment_start_indices=normalstrain_extrema_indices[:-1]
    normalstrain_segment_end_indices=normalstrain_extrema_indices[1:]
    
    normalstrain_segment_start_strain=crack_normalstrain[normalstrain_segment_start_indices]
    normalstrain_segment_end_strain=crack_normalstrain[normalstrain_segment_end_indices]

    normalstrain_maxPP = np.max(np.abs(normalstrain_segment_start_strain-normalstrain_segment_end_strain))

    normalstrain_segment_start_time=trange[normalstrain_segment_start_indices]
    normalstrain_segment_end_time=trange[normalstrain_segment_end_indices]

    normalstrain_segment_timeperiod=normalstrain_segment_end_time - normalstrain_segment_start_time

    normalstrain_segment_frequency=1.0/(2.0*normalstrain_segment_timeperiod) # A segment corresponds to peak->valley or half of a 360deg. wave
        
    
    
    # Train 1D surrogate for crack heating
    surrogate_npoints=30
    normalstrain_PP_range=np.linspace(0,normalstrain_maxPP,surrogate_npoints)
    
    normalheating_power_per_m2_hertz_side1 = np.zeros((normalstrain_PP_range.shape[0],xrange.shape[0]),dtype='d')
    normalheating_power_per_m2_hertz_stddev_side1 = np.zeros((normalstrain_PP_range.shape[0],xrange.shape[0]),dtype='d')
    
    normalheating_power_per_m2_hertz_side2 = np.zeros((normalstrain_PP_range.shape[0],xrange.shape[0]),dtype='d')
    normalheating_power_per_m2_hertz_stddev_side2 = np.zeros((normalstrain_PP_range.shape[0],xrange.shape[0]),dtype='d')
    
    for strainidx in range(normalstrain_PP_range.shape[0]):
        # We'll assume here that the static_load shift due to the mean of upper and lower)
        # is irrelevant... Really it ought to be another surrogate parameter...
    
        normal_stress_PP = plane_stress_modulus * normalstrain_PP_range[strainidx] # peak-to-peak stress
        
        (normalheating_power_per_m2_hertz_side1[strainidx,:],
         normalheating_power_per_m2_hertz_stddev_side1[strainidx,:],
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
                                                       normal_stress_PP/2.0,#vib_normal_stress_ampl.... divide by 2 converts from Peak-to-peak to amplitude
                                                       0.0, #vib_shear_stress_ampl,
                                                       1.0, #dc_excitation_frequency_numericunits.value("Hz")... we put in 1.0 for frequency and get power per Hz
                                                       crack_model_normal,
                                                       crack_model_shear,
                                                       1.0,
                                                       msqrtR,
                                                       verbose,
                                                       doplots)

        (normalheating_power_per_m2_hertz_side2[strainidx,:],
         normalheating_power_per_m2_hertz_stddev_side2[strainidx,:],
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
                                                       normal_stress_PP/2.0, #vib_normal_stress_ampl,
                                                       0.0, #vib_shear_stress_ampl,
                                                       1.0, #dc_excitation_frequency_numericunits.value("Hz")... we put in 1.0 for frequency and get power per Hz
                                                       crack_model_normal,
                                                       crack_model_shear,
                                                       1.0,
                                                       msqrtR,
                                                       verbose,
                                                       doplots)
        pass
         
    # !!!*** Need to consider uncertainty due to uncertainty in mu, msqrtR, plus response variability in here!!!***
         
    normalheating_power_per_m2_hertz_side1_surr_tck = fit_surrogate(normalstrain_PP_range,normalheating_power_per_m2_hertz_side1)
    normalheating_power_per_m2_hertz_stddev_side1_surr_tck = fit_surrogate(normalstrain_PP_range,normalheating_power_per_m2_hertz_stddev_side1)
         
    normalheating_power_per_m2_hertz_side2_surr_tck = fit_surrogate(normalstrain_PP_range,normalheating_power_per_m2_hertz_side2)
    normalheating_power_per_m2_hertz_stddev_side2_surr_tck = fit_surrogate(normalstrain_PP_range,normalheating_power_per_m2_hertz_stddev_side2)
         
    #(totalpower_side1, totalpower_stddev_side1) = integrate_power(xrange,power_per_m2_side1,power_per_m2_stddev_side1)
    #(totalpower_side2, totalpower_stddev_side2) = integrate_power(xrange,power_per_m2_side2,power_per_m2_stddev_side2)
         
    # Now calculate heating for each normalstrain_segment (via
    # normalstrain_segment_start_strain and normalstrain_segment_end_strain), iterating over x position.

    normalheating_segment_power_per_m2_side1=np.zeros((normalstrain_num_extrema,xrange.shape[0]),dtype='d')
    normalheating_segment_power_per_m2_stddev_side1=np.zeros((normalstrain_num_extrema,xrange.shape[0]),dtype='d')
    normalheating_segment_power_per_m2_side2=np.zeros((normalstrain_num_extrema,xrange.shape[0]),dtype='d')
    normalheating_segment_power_per_m2_stddev_side2=np.zeros((normalstrain_num_extrema,xrange.shape[0]),dtype='d')
    for xidx in range(xrange.shape[0]):
         # vectorized over time-segments
         straindiff = normalstrain_segment_end_strain - normalstrain_segment_start_strain
         
         normalheating_segment_power_per_m2_hertz_side1 = splev(straindiff,normalheating_power_per_m2_hertz_side1_surr_tck[xidx])
         normalheating_segment_power_per_m2_hertz_stddev_side1 = splev(straindiff,normalheating_power_per_m2_hertz_stddev_side1_surr_tck[xidx])

         normalheating_segment_power_per_m2_hertz_side2 = splev(straindiff,normalheating_power_per_m2_hertz_side2_surr_tck[xidx])
         normalheating_segment_power_per_m2_hertz_stddev_side2 = splev(straindiff,normalheating_power_per_m2_hertz_stddev_side2_surr_tck[xidx])
         
         normalheating_segment_power_per_m2_side1[:,xidx] = normalheating_segment_power_per_m2_hertz_side1*normalstrain_segment_frequency
         normalheating_segment_power_per_m2_stddev_side1[:,xidx] = normalheating_segment_power_per_m2_hertz_stddev_side1*normalstrain_segment_frequency

         normalheating_segment_power_per_m2_side2[:,xidx] = normalheating_segment_power_per_m2_hertz_side2*normalstrain_segment_frequency
         normalheating_segment_power_per_m2_stddev_side2[:,xidx] = normalheating_segment_power_per_m2_hertz_stddev_side2*normalstrain_segment_frequency
         
         pass
    
         
    normalheatingtable_power_per_m2_fine_side1 = np.zeros((trange.shape[0],xrange.shape[0]),dtype='d')
    normalheatingtable_power_per_m2_stddev_fine_side1 = np.zeros((trange.shape[0],xrange.shape[0]),dtype='d')
    normalheatingtable_power_per_m2_fine_side2 = np.zeros((trange.shape[0],xrange.shape[0]),dtype='d')
    normalheatingtable_power_per_m2_stddev_fine_side2 = np.zeros((trange.shape[0],xrange.shape[0]),dtype='d')
    
    for timesegidx in range(normalstrain_num_extrema):
         normalheatingtable_power_per_m2_fine_side1[normalstrain_segment_start_indices[timesegidx]:normalstrain_segment_end_indices[timesegidx],:]=normalheating_segment_power_per_m2_side1[timesegidx,:]
         normalheatingtable_power_per_m2_stddev_fine_side1[normalstrain_segment_start_indices[timesegidx]:normalstrain_segment_end_indices[timesegidx],:]=normalheating_segment_power_per_m2_stddev_side1[timesegidx,:]
         normalheatingtable_power_per_m2_fine_side2[normalstrain_segment_start_indices[timesegidx]:normalstrain_segment_end_indices[timesegidx],:]=normalheating_segment_power_per_m2_side2[timesegidx,:]
         normalheatingtable_power_per_m2_stddev_fine_side2[normalstrain_segment_start_indices[timesegidx]:normalstrain_segment_end_indices[timesegidx],:]=normalheating_segment_power_per_m2_stddev_side2[timesegidx,:]
         pass

    
    # Need to filter down the fine time data to coarser time grid
    normalheatingtable_power_per_m2_side1 = filterfine(normalheatingtable_power_per_m2_fine_side1,trange,t)
    normalheatingtable_power_per_m2_side2 = filterfine(normalheatingtable_power_per_m2_fine_side2,trange,t)

    # normalstrain_heating = surrogate_calc_heating()
    
    meanpower_per_m2_side1 = np.mean(normalheatingtable_power_per_m2_side1,axis=0) # average over time... should add in shear
    meanpower_per_m2_side2 = np.mean(normalheatingtable_power_per_m2_side2,axis=0) # average over time... should add in shear
         
    totalpower_side1 = integrate_power(xrange,meanpower_per_m2_side1)
    totalpower_side2 = integrate_power(xrange,meanpower_per_m2_side2)
    

    totalpower=totalpower_side1 + totalpower_side2
     
    pl.figure()
    pl.clf()
    pl.plot(-xrange*1e3,meanpower_per_m2_side1/1.e3,'-',
            xrange*1e3,meanpower_per_m2_side2/1.e3,'-',)
    pl.grid()
    pl.xlabel('Radius from center (mm)')
    pl.ylabel('Heating power (kW/m^2)')
    pl.title('Average crack power deposition')

    heatpower_href = hrefv(quote(dc_measident_str+"_heatpower.png"),dc_dest_href)
    pl.savefig(heatpower_href.getpath(),dpi=300)

    heatingfh = open(dc_heatingdata_href.getpath(),"w")
    heatingfh.write("% t(s) \t r(m) \t side1_heating(W/m^2) \t side2_heating(W/m^2)\n")

    for tidx in range(t.shape[0]):
        for xidx in range(xrange.shape[0]):
            heatingfh.write("%.8e\t%.8e\t%.8e\t%.8e\n" % (t[tidx],xrange[xidx],normalheatingtable_power_per_m2_side1[tidx,xidx],normalheatingtable_power_per_m2_side2[tidx,xidx]))
            pass
        pass
    heatingfh.close()
            
    

    ret= { 
        "dc:closureplot_side1": closure_state_side1_href,
        "dc:closureplot_side2": closure_state_side2_href,
        "dc:heatpower": heatpower_href,  # heating plot
        "dc:heatingdata": dc_heatingdata_href, # TSV (tab separated value) data for COMSOL
        "dc:heatingtotalpower": numericunitsv(totalpower,"W"),
    }

    return ret
