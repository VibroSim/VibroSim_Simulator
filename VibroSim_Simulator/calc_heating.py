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


from crackclosuresim2 import inverse_closure
from crackclosuresim2 import crackopening_from_tensile_closure
from crackclosuresim2 import solve_normalstress

from crackclosuresim2 import ModeI_throughcrack_CODformula
from crackclosuresim2 import Tada_ModeI_CircularCrack_along_midline
from crackclosuresim2 import ModeII_throughcrack_CSDformula
from crackclosuresim2.fabrikant import Fabrikant_ModeII_CircularCrack_along_midline

from angled_friction_model.angled_friction_model import angled_friction_model
from angled_friction_model.angled_friction_model import integrate_power




def fit_surrogate(strainvals,heatingvals):
    """ Fit a 1D surrogate over axis 0 of 2D array heatingvals, Each row of heatingvals corresponsd to
    the corresponding entry in strainvals. The operation is vectorized over the columns (axis 1) of 
    heatingvals"""
    
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
    """Use a running average filter along axis 0 (down the rows) of finedata,
    which is sampled according to the timebase given in finetrange.

    Then downsample and return points sampled according to the timebase
    given in coarsetrange.

    The number of points in the running average filter is roughly 
    2*coarse_timestep/fine_timestep
"""
    
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


def strain_segments_between_extrema(trange,crack_strain):
    """ Break crack_strain, a function of time with timebase trange, 
    into segments connecting extrema -- the first point, all the relative 
    extrema in the waveform, and the last point. Return the characteristics
    of each segment"""

    # Find extrema
    strain_relative_maxima = (crack_strain[1:-1] > crack_strain[2:]) & (crack_strain[1:-1] > crack_strain[:-2])
    strain_relative_minima = (crack_strain[1:-1] < crack_strain[2:]) & (crack_strain[1:-1] < crack_strain[:-2])
    strain_extrema = strain_relative_minima | strain_relative_maxima
    strain_extrema_indices = np.concatenate(((0,),np.where(strain_extrema)[0]+1,(crack_strain.shape[0]-1,)))
    
    strain_num_extrema=strain_extrema_indices.shape[0]-1

    strain_segment_start_indices=strain_extrema_indices[:-1]
    strain_segment_end_indices=strain_extrema_indices[1:]
    
    strain_segment_start_strain=crack_strain[strain_segment_start_indices]
    strain_segment_end_strain=crack_strain[strain_segment_end_indices]

    strain_maxPP = np.max(np.abs(strain_segment_start_strain-strain_segment_end_strain))

    strain_segment_start_time=trange[strain_segment_start_indices]
    strain_segment_end_time=trange[strain_segment_end_indices]

    strain_segment_timeperiod=strain_segment_end_time - strain_segment_start_time

    strain_segment_frequency=1.0/(2.0*strain_segment_timeperiod) # A segment corresponds to peak->valley or half of a 360deg. wave

    return (strain_maxPP,  # Maximum peak-to-peak variation in strain
            strain_segment_start_indices,
            strain_segment_end_indices,
            strain_segment_start_strain,
            strain_segment_end_strain,
            strain_segment_frequency)



def train_surrogate(surrogate_type, # either "normal" or "shear"
                    surrogate_npoints,
                    strain_maxPP,
                    x_bnd,xrange,xstep, # range of radii from crack center
                    numdraws,
                    YoungsModulus, # Pa
                    PoissonsRatio,
                    sigma_yield,
                    friction_coefficient,
                    msqrtR,
                    closure_stress, # closure stress field for this side of crack
                    crack_opening, # crack opening displacement for this side of crack
                    angular_stddev, # tortuosity, in radians
                    aside, # length of this side of crack,
                    static_load, # static load (if any) to be superimposed on preexisting closure stress and crack opening
                    crack_model_normal_name,
                    crack_model_shear_name,
                    crack_model_shear_factor): # shear sensitivity factor (nominally 1.0)
    
    verbose=False
    doplots=False

    
    
    plane_stress_modulus = YoungsModulus/(1.0-PoissonsRatio**2.0)
    tau_yield = sigma_yield/2.0
    
    crack_model_normal = crack_model_normal_by_name(crack_model_normal_name,YoungsModulus,PoissonsRatio)
    crack_model_shear = crack_model_shear_by_name(crack_model_shear_name,YoungsModulus,PoissonsRatio)
    
    


    # Train 1D surrogate for crack heating
    strain_PP_range=np.linspace(0,strain_maxPP,surrogate_npoints)
    
    heating_power_per_m2_hertz = np.zeros((strain_PP_range.shape[0],xrange.shape[0]),dtype='d')
    heating_power_per_m2_hertz_stddev = np.zeros((strain_PP_range.shape[0],xrange.shape[0]),dtype='d')
    
    for strainidx in range(strain_PP_range.shape[0]):
        # We'll assume here that the static_load shift due to the mean of upper and lower)
        # is irrelevant... Really it ought to be another surrogate parameter...
    
        stress_PP = plane_stress_modulus * strain_PP_range[strainidx] # peak-to-peak stress

        if surrogate_type == "normal":
            vib_normal_stress_ampl = stress_PP/2.0 #... divide by 2 converts from Peak-to-peak to amplitude
            vib_shear_stress_ampl = 0.0
            pass
        elif surrogate_type == "shear":
            vib_normal_stress_ampl = 0.0
            vib_shear_stress_ampl = stress_PP/2.0 #... divide by 2 converts from Peak-to-peak to amplitude
            pass
        else:
            raise ValueError("Unknown surrogate type: %s" % (surrogate_type))
        
        (heating_power_per_m2_hertz[strainidx,:],
         heating_power_per_m2_hertz_stddev[strainidx,:],
         vibration_ampl) = angled_friction_model(x_bnd,xrange,xstep,
                                                       numdraws,
                                                       YoungsModulus,
                                                       PoissonsRatio,
                                                       sigma_yield,tau_yield,
                                                       friction_coefficient,
                                                       closure_stress,
                                                       crack_initial_opening,
                                                       angular_stddev,
                                                       aside,
                                                       static_load,
                                                       vib_normal_stress_ampl,
                                                       vib_shear_stress_ampl,
                                                       1.0, #dc_excitation_frequency_numericunits.value("Hz")... we put in 1.0 for frequency and get power per Hz
                                                       crack_model_normal,
                                                       crack_model_shear,
                                                       crack_model_shear_factor,
                                                       msqrtR,
                                                       verbose,
                                                       doplots)

        pass
         
    heating_power_per_m2_hertz_surr_tck = fit_surrogate(strain_PP_range,heating_power_per_m2_hertz)
    heating_power_per_m2_hertz_stddev_surr_tck = fit_surrogate(strain_PP_range,heating_power_per_m2_hertz_stddev)
         
    return (heating_power_per_m2_hertz_surr_tck,heating_power_per_m2_hertz_stddev_surr_tck)


def evaluate_per_segment_heating_power(xrange,
                                       heating_power_per_m2_hertz_surr_tck, # array of (tck) tuples with surrogates for heating power, one surrogate per x position
                                       heating_power_per_m2_hertz_stddev_surr_tck, # array of (tck) tuples with surrogates for standard deviation, one surrogate per x position
                                       segment_start_strain,
                                       segment_end_strain,
                                       segment_frequency):
    
    num_extrema = segment_start_strain.shape[0]
    
    heating_segment_power_per_m2=np.zeros((num_extrema,xrange.shape[0]),dtype='d')
    heating_segment_power_per_m2_stddev=np.zeros((num_extrema,xrange.shape[0]),dtype='d')

    straindiff = np.abs(segment_end_strain - segment_start_strain)
    
    for xidx in range(xrange.shape[0]):
         # vectorized over time-segments
         
        heating_segment_power_per_m2_hertz = splev(straindiff,heating_power_per_m2_hertz_surr_tck[xidx],ext=2)
        heating_segment_power_per_m2_hertz[heating_segment_power_per_m2_hertz < 0.0] = 0.0 # negative heating not allowed (in case spline fit goes negative for some reason)
        heating_segment_power_per_m2_hertz_stddev = splev(straindiff,heating_power_per_m2_hertz_stddev_surr_tck[xidx],ext=2)
        heating_segment_power_per_m2_hertz_stddev[heating_segment_power_per_m2_hertz_stddev < 0.0] = 0.0 # negative stddev not allowed (in case spline fit goes negative for some reason)
        
        heating_segment_power_per_m2[:,xidx] = heating_segment_power_per_m2_hertz*segment_frequency
        heating_segment_power_per_m2_stddev[:,xidx] = heating_segment_power_per_m2_hertz_stddev*segment_frequency
        
        pass
    return (
        heating_segment_power_per_m2,
        heating_segment_power_per_m2_stddev
    )

def heatingpower_timeseries_from_segments(trange,xrange,segment_start_indices,segment_end_indices,segment_values):
    """
    segment_values has shape(num_segments,xrange.shape[0])
    segment_values_fine will have shape(trange.shape[0],xrange.shape[0])
    segment_start_indices and segment_end_indices define the starting and ending index into trange for each segment. 
    
    This routine fills the segment values into the resulting segment_values_fine based on the start and end indices"""

    segment_values_fine = np.zeros((trange.shape[0],xrange.shape[0]),dtype='d')

    num_extrema = segment_start_indices.shape[0]
    
    for timesegidx in range(num_extrema):
        segment_values_fine[segment_start_indices[timesegidx]:segment_end_indices[timesegidx],:]=segment_values[timesegidx,:]
        pass

    return segment_values_fine

def calc_heating_welder(friction_coefficient,
                        msqrtR,  # m^1.5
                        static_load, # Pa
                        angular_stddev, # Tortuosity in radians
                        numdraws, 
                        YoungsModulus, # Pa
                        PoissonsRatio, 
                        sigma_yield, # Pa
                        Density, # kg/m^3                        
                        x_bnd,xrange,xstep, # base of crack radius values
                        closure_stress_side1, # side1 closure state
                        crack_opening_side1,
                        aside1,                        
                        closure_stress_side2, # side2 clsoure state
                        crack_opening_side2,
                        aside2,                        
                        crack_model_normal_name,
                        crack_model_shear_name,
                        crack_model_shear_factor,
                        exc_t0,exc_t4, # s
                        motion_path,
                        heatingdata_path):


    # note t and tstep are the steps for heating output table, 
    # as distinct from trange and dt which correspond to the motion file
    # and are much finer
    t=np.linspace(exc_t0,exc_t4,100)
    t_step = t[1]-t[0]
    
    
    shear_modulus = YoungsModulus/(2.0*(1.0+PoissonsRatio))

    # Load motion file
    motiontable = pd.read_csv(motion_path,index_col=0)
                        
                        
    # Note distinction between trange and dt defined here and the 
    # t and t_step defined for the heating output table above. 
    trange = np.array(motiontable.index) # ["Time(s)"]
    dt=trange[1]-trange[0]
                        

    # ***!!! Should crosscheck with single frequency result by evaluating a motion table from single frequency data
    # and running that through this !!!***
    
    crack_normalstrain = np.array(motiontable["specimen_crackcenternormalstrain"])
    crack_shearstrain = np.array(motiontable["specimen_crackcentershearstrain"])

    crack_strain_fig = pl.figure()
    pl.plot(trange*1e3,crack_normalstrain*1e6,'-',
            trange*1e3,crack_shearstrain*1e6,'-')
    pl.grid(True)
    pl.legend(('Normal strain','shear strain'))
    pl.xlabel('Time (ms)')
    pl.ylabel('Engineering strain (micros)')
    pl.title('Crack strain')


    
    
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





    (normalstrain_maxPP,  # Maximum peak-to-peak variation in strain
     normalstrain_segment_start_indices,
     normalstrain_segment_end_indices,
     normalstrain_segment_start_strain,
     normalstrain_segment_end_strain,
     normalstrain_segment_frequency) = strain_segments_between_extrema(trange,crack_normalstrain)
    
    (shearstrain_maxPP,  # Maximum peak-to-peak variation in strain
     shearstrain_segment_start_indices,
     shearstrain_segment_end_indices,
     shearstrain_segment_start_strain,
     shearstrain_segment_end_strain,
     shearstrain_segment_frequency) = strain_segments_between_extrema(trange,crack_shearstrain)
                        
    

                        
    


    # !!!*** Need to consider uncertainty due to uncertainty in mu, msqrtR, plus response variability in here!!!***
         
    (
        normalheating_power_per_m2_hertz_side1_surr_tck,
        normalheating_power_per_m2_hertz_stddev_side1_surr_tck
    ) = train_surrogate("normal",
                        30, # surrogate_npoints
                        normalstrain_maxPP,
                        x_bnd,xrange,xstep,
                        numdraws
                        YoungsModulus, # Pa
                        PoissonsRatio,
                        sigma_yield,
                        friction_coefficient,
                        msqrtR,
                        closure_stress_side1,
                        crack_opening_side1,
                        angular_stddev, # tortuosity, in radians,
                        aside1, # length of side of crack
                        static_load,
                        crack_model_normal_name,
                        crack_model_shear_name,
                        crack_model_shear_factor) # shear sensitivity factor (nominally 1.0)

    
    (
        normalheating_power_per_m2_hertz_side2_surr_tck,
        normalheating_power_per_m2_hertz_stddev_side2_surr_tck
    ) = train_surrogate("normal",
                        30, # surrogate_npoints
                        normalstrain_maxPP,
                        x_bnd,xrange,xstep,
                        numdraws
                        YoungsModulus, # Pa
                        PoissonsRatio,
                        sigma_yield,
                        friction_coefficient,
                        msqrtR,
                        closure_stress_side2,
                        crack_opening_side2,
                        angular_stddev, # tortuosity, in radians,
                        aside2, # length of side of crack
                        static_load,
                        crack_model_normal_name,
                        crack_model_shear_name,
                        crack_model_shear_factor) # shear sensitivity factor (nominally 1.0)


    
    (
        shearheating_power_per_m2_hertz_side1_surr_tck,
        shearheating_power_per_m2_hertz_stddev_side1_surr_tck
    ) = train_surrogate("shear",
                        30, # surrogate_npoints
                        shearstrain_maxPP,
                        x_bnd,xrange,xstep,
                        numdraws
                        YoungsModulus, # Pa
                        PoissonsRatio,
                        sigma_yield,
                        friction_coefficient,
                        msqrtR,
                        closure_stress_side1,
                        crack_opening_side1,
                        angular_stddev, # tortuosity, in radians,
                        aside1, # length of side of crack
                        static_load,
                        crack_model_normal_name,
                        crack_model_shear_name,
                        crack_model_shear_factor) # shear sensitivity factor (nominally 1.0)
    
    (
        shearheating_power_per_m2_hertz_side2_surr_tck,
        shearheating_power_per_m2_hertz_stddev_side2_surr_tck
    ) = train_surrogate("shear",
                        30, # surrogate_npoints
                        shearstrain_maxPP,
                        x_bnd,xrange,xstep,
                        numdraws
                        YoungsModulus, # Pa
                        PoissonsRatio,
                        sigma_yield,
                        friction_coefficient,
                        msqrtR,
                        closure_stress_side2,
                        crack_opening_side2,
                        angular_stddev, # tortuosity, in radians,
                        aside2, # length of side of crack
                        static_load,
                        crack_model_normal_name,
                        crack_model_shear_name,
                        crack_model_shear_factor) # shear sensitivity factor (nominally 1.0)

         
    #(totalpower_side1, totalpower_stddev_side1) = integrate_power(xrange,power_per_m2_side1,power_per_m2_stddev_side1)
    #(totalpower_side2, totalpower_stddev_side2) = integrate_power(xrange,power_per_m2_side2,power_per_m2_stddev_side2)
         
    # Now calculate heating for each normalstrain_segment (via
    # normalstrain_segment_start_strain and normalstrain_segment_end_strain), iterating over x position.

    (
        normalheating_segment_power_per_m2_side1,
        normalheating_segment_power_per_m2_stddev_side1
    ) = evaluate_per_segment_heating_power(xrange,
                                           normalheating_power_per_m2_hertz_side1_surr_tck, # array of (tck) tuples with surrogates for heating power one surrogate per x position
                                           normalheating_power_per_m2_hertz_stddev_side1_surr_tck, # array of (tck) tuples with surrogates for standard deviation, one surrogate per x position
                                           normalstrain_segment_start_strain,
                                           normalstrain_segment_end_strain,
                                           normalstrain_segment_frequency)
    
    (
        normalheating_segment_power_per_m2_side2,
        normalheating_segment_power_per_m2_stddev_side2
    ) = evaluate_per_segment_heating_power(xrange,
                                           normalheating_power_per_m2_hertz_side2_surr_tck, # array of (tck) tuples with surrogates for heating power one surrogate per x position
                                           normalheating_power_per_m2_hertz_stddev_side2_surr_tck, # array of (tck) tuples with surrogates for standard deviation, one surrogate per x position
                                           normalstrain_segment_start_strain,
                                           normalstrain_segment_end_strain,
                                           normalstrain_segment_frequency)


    (
        shearheating_segment_power_per_m2_side1,
        shearheating_segment_power_per_m2_stddev_side1
    ) = evaluate_per_segment_heating_power(xrange,
                                           shearheating_power_per_m2_hertz_side1_surr_tck, # array of (tck) tuples with surrogates for heating power one surrogate per x position
                                           shearheating_power_per_m2_hertz_stddev_side1_surr_tck, # array of (tck) tuples with surrogates for standard deviation, one surrogate per x position
                                           shearstrain_segment_start_strain,
                                           shearstrain_segment_end_strain,
                                           shearstrain_segment_frequency)
    
    (
        shearheating_segment_power_per_m2_side2,
        shearheating_segment_power_per_m2_stddev_side2
    ) = evaluate_per_segment_heating_power(xrange,
                                           shearheating_power_per_m2_hertz_side2_surr_tck, # array of (tck) tuples with surrogates for heating power one surrogate per x position
                                           shearheating_power_per_m2_hertz_stddev_side2_surr_tck, # array of (tck) tuples with surrogates for standard deviation, one surrogate per x position
                                           shearstrain_segment_start_strain,
                                           shearstrain_segment_end_strain,
                                           shearstrain_segment_frequency)
    

    normalheatingtable_power_per_m2_fine_side1 = heatingpower_timeseries_from_segments(trange,xrange,
                                                                                       normalstrain_segment_start_indices,
                                                                                       normalstrain_segment_end_indices,
                                                                                       normalheating_segment_power_per_m2_side1)
    normalheatingtable_power_per_m2_stddev_fine_side1 = heatingpower_timeseries_from_segments(trange,xrange,
                                                                                              normalstrain_segment_start_indices,
                                                                                              normalstrain_segment_end_indices,
                                                                                              normalheating_segment_power_per_m2_stddev_side1)
    
    
    normalheatingtable_power_per_m2_fine_side2 = heatingpower_timeseries_from_segments(trange,xrange,
                                                                                       normalstrain_segment_start_indices,
                                                                                       normalstrain_segment_end_indices,
                                                                                       normalheating_segment_power_per_m2_side2)
    
    normalheatingtable_power_per_m2_stddev_fine_side2 = heatingpower_timeseries_from_segments(trange,xrange,
                                                                                              normalstrain_segment_start_indices,
                                                                                              normalstrain_segment_end_indices,
                                                                                              normalheating_segment_power_per_m2_stddev_side2)

    
    shearheatingtable_power_per_m2_fine_side1 = heatingpower_timeseries_from_segments(trange,xrange,
                                                                                      shearstrain_segment_start_indices,
                                                                                      shearstrain_segment_end_indices,
                                                                                      shearheating_segment_power_per_m2_side1)
    
    shearheatingtable_power_per_m2_stddev_fine_side1 = heatingpower_timeseries_from_segments(trange,xrange,
                                                                                             shearstrain_segment_start_indices,
                                                                                             shearstrain_segment_end_indices,
                                                                                             shearheating_segment_power_per_m2_stddev_side1)

    
    shearheatingtable_power_per_m2_fine_side2 = heatingpower_timeseries_from_segments(trange,xrange,
                                                                                      shearstrain_segment_start_indices,
                                                                                      shearstrain_segment_end_indices,
                                                                                      shearheating_segment_power_per_m2_side2)

    shearheatingtable_power_per_m2_stddev_fine_side2 = heatingpower_timeseries_from_segments(trange,xrange,
                                                                                             shearstrain_segment_start_indices,
                                                                                             shearstrain_segment_end_indices,
                                                                                             shearheating_segment_power_per_m2_stddev_side2)

    
    
    # Need to filter down the fine time data to coarser time grid
    normalheatingtable_power_per_m2_side1 = filterfine(normalheatingtable_power_per_m2_fine_side1,trange,t)
    normalheatingtable_power_per_m2_side2 = filterfine(normalheatingtable_power_per_m2_fine_side2,trange,t)

    shearheatingtable_power_per_m2_side1 = filterfine(shearheatingtable_power_per_m2_fine_side1,trange,t)
    shearheatingtable_power_per_m2_side2 = filterfine(shearheatingtable_power_per_m2_fine_side2,trange,t)

    # (Should provide corresponding stddev estimates) 



    normal_heatgram_side1_fig=pl.figure()
    pl.clf()
    pl.imshow(normalheatingtable_power_per_m2_side1,aspect='auto',origin='lower',extent=((xrange[0]-xstep/2.0)*1e3,(xrange[-1]+xstep/2.0)*1e3,t[0]-t_step/2.0,t[-1]+t_step/2.0))
    pl.colorbar()
    pl.grid(True)
    pl.xlabel('Radius from crack center (mm)')
    pl.ylabel('Time (s)')
    pl.title('Heating power due to normal strain (W/m^2), side1')


    normal_heatgram_side2_fig=pl.figure()
    pl.clf()
    pl.imshow(normalheatingtable_power_per_m2_side2,aspect='auto',origin='lower',extent=((xrange[0]-xstep/2.0)*1e3,(xrange[-1]+xstep/2.0)*1e3,t[0]-t_step/2.0,t[-1]+t_step/2.0))
    pl.colorbar()
    pl.grid(True)
    pl.xlabel('Radius from crack center (mm)')
    pl.ylabel('Time (s)')
    pl.title('Heating power due to normal strain (W/m^2), side2')
    


        
    meanpower_per_m2_side1 = np.mean(normalheatingtable_power_per_m2_side1,axis=0) + np.mean(shearheatingtable_power_per_m2_side1,axis=0) # average over time... 
    meanpower_per_m2_side2 = np.mean(normalheatingtable_power_per_m2_side2,axis=0) + np.mean(shearheatingtable_power_per_m2_side2,axis=0) # average over time...
         
    totalpower_side1 = integrate_power(xrange,meanpower_per_m2_side1)
    totalpower_side2 = integrate_power(xrange,meanpower_per_m2_side2)
    

    totalpower=totalpower_side1 + totalpower_side2
     

    heatpower_fig=pl.figure()
    pl.clf()
    pl.plot(-xrange*1e3,meanpower_per_m2_side1/1.e3,'-',
            xrange*1e3,meanpower_per_m2_side2/1.e3,'-',)
    pl.grid()
    pl.xlabel('Radius from center (mm)')
    pl.ylabel('Heating power (kW/m^2)')
    pl.title('Average crack power deposition')


    heatingfh = open(heatingdata_path,"w")
    heatingfh.write("% t(s) \t r(m) \t side1_heating(W/m^2) \t side2_heating(W/m^2)\n")
    
    for tidx in range(t.shape[0]):
        for xidx in range(xrange.shape[0]):
            heatingfh.write("%.8e\t%.8e\t%.8e\t%.8e\n" % (t[tidx],xrange[xidx],normalheatingtable_power_per_m2_side1[tidx,xidx] + shearheatingtable_power_per_m2_side1[tidx,xidx],normalheatingtable_power_per_m2_side2[tidx,xidx] + shearheatingtable_power_per_m2_side2[tidx,xidx]))
            pass
        pass
    heatingfh.close()
    

    
    return (crack_strain_fig,
            normal_heatgram_side1_fig,
            normal_heatgram_side2_fig,
            shear_heatgram_side1_fig,
            shear_heatgram_side2_fig,
            heatpower_fig,
            normalheatingtable_power_per_m2_side1,
            normalheatingtable_power_per_m2_side2,
            shearheatingtable_power_per_m2_side1,
            shearheatingtable_power_per_m2_side2,
            meanpower_per_m2_side1
            meanpower_per_m2_side2
            totalpower)
            
            



def calc_heating_singlefrequency(friction_coefficient,
                                 msqrtR,  # m^1.5
                                 static_load, # Pa
                                 angular_stddev, # Tortuosity in radians
                                 numdraws, 
                                 YoungsModulus, # Pa
                                 PoissonsRatio, 
                                 sigma_yield, # Pa
                                 Density, # kg/m^3                        
                                 x_bnd,xrange,xstep, # base of crack radius values
                                 closure_stress_side1, # side1 closure state
                                 crack_opening_side1,
                                 aside1,                        
                                 closure_stress_side2, # side2 clsoure state
                                 crack_opening_side2,
                                 aside2,                        
                                 crack_model_normal_name,
                                 crack_model_shear_name,
                                 crack_model_shear_factor,
                                 exc_t0,exc_t3, # s
                                 harmonicburst_normalstrain, # NOTE: complex
                                 harmonicburst_shearstrain, # NOTE: complex                
                                 heatingdata_path):


    # note t and tstep are the steps for heating output table, 
    # as distinct from trange and dt which correspond to the motion file
    # and are much finer
    t=np.linspace(exc_t0,exc_t4,100)
    t_step = t[1]-t[0]
    
    
    plane_stress_modulus = YoungsModulus/(1.0-PoissonsRatio**2.0)
    shear_modulus = YoungsModulus/(2.0*(1.0+PoissonsRatio))
    tau_yield = sigma_yield/2.0
    
                        
    # Note distinction between trange and dt defined here and the 
    # t and t_step defined for the heating output table above. 
    trange = np.array(motiontable.index) # ["Time(s)"]
    dt=trange[1]-trange[0]
                        
    
    # Note: We are ignoring the phases of the normal and shear strains
    vib_normal_stress_ampl = plane_stress_modulus * abs(harmonicburst_normalstrain)
    vib_shear_stress_ampl = shear_modulus * abs(harmonicburst_shearstrain)

    crack_model_normal = crack_model_normal_by_name(crack_model_normal_name,YoungsModulus,PoissonsRatio)
    crack_model_shear = crack_model_shear_by_name(crack_model_shear_name,YoungsModulus,PoissonsRatio)

    
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


    
    # !!!*** Need to consider uncertainty due to uncertainty in mu, msqrtR, plus response variability in here!!!***

    totalpower=totalpower_side1 + totalpower_side2




    heatpower_fig = pl.figure()
    pl.clf()
    pl.plot(-xrange*1e3,power_per_m2_side1/1.e3,'-',
            xrange*1e3,power_per_m2_side2/1.e3,'-',)
    pl.grid()
    pl.xlabel('Radius from center (mm)')
    pl.ylabel('Heating power (kW/m^2)')
    pl.title('Crack power deposition')


    heatingfh = open(heatingdata_path,"w")
    heatingfh.write("% t(s) \t r(m) \t side1_heating(W/m^2) \t side2_heating(W/m^2)\n")
    
    for tidx in range(t.shape[0]):
        for xidx in range(xrange.shape[0]):
            heatingfh.write("%.8e\t%.8e\t%.8e\t%.8e\n" % (t[tidx],xrange[xidx],power_per_m2_side1[xidx],power_per_m2_side2[xidx]))
            pass
        pass
    heatingfh.close()
            

    return (heatpower_fig,
            power_per_m2_side1,
            power_per_m2_side2,
            totalpower)
