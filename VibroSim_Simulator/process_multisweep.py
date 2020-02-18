
import os
import os.path
import sys
import numpy as np
import scipy.interpolate
import collections
from . import read_comsol_probe_txt
#import lzma
import pandas as pd
#import bz2

from matplotlib import pyplot as pl

def array2string(array):
    """... because nparray2string is way too slow for million+ element arrays
    ...similar to np.array2string(array,separator=',\n',suppress_small=False,threshold=np.inf,floatmode='unique')"""
    
    el_reprs = [ repr(el) for el in array ]
    return ",\n".join(el_reprs)



def interpolate_spectrum(min_df,dt,FreqResp,impulseexcitation_width):
    
    
    max_t = 1.0/min_df
    n = int(round(max_t/dt))
    if False and not((n & (n-1) == 0) and n != 0):
        # n is not a power of 2... we don't rely on this here
        # but the way we generated the steps (currently) always
        # generates a power of 2. So raise an error flag because
        # this suggests a problem! -- perhaps our assumed dt doesn't match
        raise ValueError("Number of steps should be a power of 2!")


    InputFreqs = np.array(list(FreqResp.keys()),dtype='d')
    InputMeas = np.array(list(FreqResp.values()),dtype=np.complex128)
    
    #pl.figure()
    #pl.loglog(InputFreqs,np.abs(InputVels))
    #pl.title('Magnitude')
    
    #pl.figure()
    #pl.loglog(InputFreqs,np.real(InputVels))
    #pl.title('Magnitude')

    #pl.figure()
    #pl.semilogx(InputFreqs,np.angle(InputVels))
    #pl.title('Magnitude')


    
    trange = np.arange(n,dtype='d')*dt
    
    df = 1.0/(n*dt)
    
    frange = np.arange(n,dtype='d')*df
    
    interpolator = scipy.interpolate.interp1d(InputFreqs,InputMeas,kind="cubic",fill_value="extrapolate")  # lowpass interpolation would be better, but scipy doesn't provide...
    
    meas=2.0*interpolator(frange)  # m/s per N or m per N; factor of 2.0 is because we are neglecting the negative frequencies (then after we IFFT we will neglect the imaginary part)
    
    meas_timedomain=np.fft.ifft(meas)*n*df # m/s per (N*s)

    #pl.figure()
    #pl.plot(trange,np.real(meas_timedomain),'-',
    #        trange,np.imag(meas_timedomain),'-')
    
    # This plot should roughly match COMSOL but is unfiltered.
    # the minus sign is because the COMSOL model calculation
    # depending on which is used, evaluates positive in the
    # outward normal direction. 
    #pl.figure()
    #pl.plot(trange,-np.real(meas_timedomain),'-')


    # Attempt convolution with a Gaussian pulse of area 1
    GPulse_sigma = impulseexcitation_width # parameter from COMSOL
    trange_shifted=trange.copy()
    trange_shifted[n//2:] -= n*dt
    GPulse = 1.0/np.sqrt(2*np.pi*GPulse_sigma**2.0)*np.exp(-trange_shifted**2.0/(2.0*GPulse_sigma**2.0))
    GPulse_fft=np.fft.fft(GPulse)*dt
    # Zero out negative frequencies -- for complex time domain signal
    GPulse_fft[n//2:]=0.0

    meas_filtered = meas*GPulse_fft
    meas_filtered_timedomain=np.fft.ifft(meas_filtered)*n*df # m/s per (N*s)
    
    return (trange,frange,meas_timedomain,meas,meas_filtered_timedomain,meas_filtered)


def read_spectrum_ansys(seg_table_filenames,seg_tables,dt,impulseexcitation_width,location,quantity,amplitude_units,amplitude_multiplier,phase_units,phase_multiplier):
    """Read spectrum from the given segment tables. 
    location, quantity, and units, joined by spaces give the column title
    to use for the data. The frequency is presumed to be in a column
    titled "Frequency (Hz)". 

    Evaluating a single specrum requires both an amplitude and a phase, each of which 
    will be multiplied by the given multiplier prior to any processing. The multiplier 
    should be the correct factor to convert from the given units to base MKS units 
    (radians for phase). """

    freq_key = "Frequency (Hz)"

    amplitude_to_join = [ component for component in (location,quantity,"Amplitude",amplitude_units) if component is not None ]
    amplitude_title = " ".join(amplitude_to_join)

    phase_to_join = [ component for component in (location,quantity,"Phase Angle",phase_units) if component is not None ]
    phase_title = " ".join(phase_to_join)

    descr_to_join = [ component for component in (location,quantity) if component is not None ]
    descr = " ".join(descr_to_join)
    
    FreqResp = collections.OrderedDict()

    dfs = []
    

    for segnum in range(len(seg_tables)):
        this_segdata = []

        if amplitude_title not in seg_tables[segnum]:
            raise ValueError("Table %s missing column \"%s\"" % (seg_table_filenames[segnum],amplitude_title))

        if phase_title not in seg_tables[segnum]:
            raise ValueError("Table %s missing column \"%s\"" % (seg_table_filenames[segnum],phase_title))

        for (index,row) in seg_tables[segnum].iterrows():
            complex_response = row[amplitude_title]*amplitude_multiplier*np.exp((0+1j)*row[phase_title]*phase_multiplier)
            FreqResp[row[freq_key]] = complex_response
            this_segdata.append(complex_response)
            pass
        
        this_segarray=np.array(this_segdata,dtype=np.complex128)

        # plot up the separate inverse transform of this segment alone, assuming all other segments are zero
        # df = 1/(n*dt)
        # dt = 1/(n*df)
        # t_total = n*dt = 1/df
        
        # interpolate to a finer dt given_dt = (1/(n_new*df)) by padding with zeros in freq domain. 

        # this_dt_uninterpolated = 1.0/(len(this_segdata)*this_df)

        this_df = seg_tables[segnum].iloc[1][freq_key] - seg_tables[segnum].iloc[0][freq_key]

        dfs.append(this_df)
        
        # given_dt = 1.0/(n_new*this_df)
        # ->  n_new = 1.0/(given_dt*this_df)
        n_new = int(np.round(1.0/(dt*this_df)))
        
        this_freqbase = seg_tables[segnum].iloc[0][freq_key] + np.arange(n_new)*this_df
        
        this_timebase = np.arange(n_new)*(1.0/(n_new*this_df))

        # multiply inverse transform by n_new (to compensate for built-in normalization of ifft) 
        # by df -- quantity integrated over -- and by exponential that compensates for the frequency 
        # corresponding to our first step
        this_segtimedata = np.fft.ifft(this_segarray,n=n_new)*n_new*this_df*np.exp(2.0*np.pi*(0+1j)*this_freqbase[0]*this_timebase)

        pl.figure()
        pl.plot(this_timebase,np.real(this_segtimedata),'-')
        pl.xlabel('Time (s)')
        pl.ylabel(descr)
        pl.title('seg%d' % (segnum+1))

        pass

    min_df = np.min(dfs)
    
    (trange,frange,meas_timedomain,meas,meas_filtered_timedomain,meas_filtered) = interpolate_spectrum(min_df,dt,FreqResp,impulseexcitation_width)
    


    return (trange,frange,meas_timedomain,meas,meas_filtered_timedomain,meas_filtered)
    


def process_multisweep_ansys(laser_name,crack_name,plotdir,plotprefix,dt,impulseexcitation_width,endcrop,seg_table_names,output_filename):
    numsegs = len(seg_table_names)

    seg_tables = [ pd.read_excel(seg_table_name) for seg_table_name in seg_table_names ]

        
    (xducer_vel_trange,xducer_velspec_frange,xducer_vel_timedomain,xducer_velspec,xducer_vel_filtered_timedomain,xducer_velspec_filtered) = read_spectrum_ansys(seg_table_names,seg_tables,dt,impulseexcitation_width,"Horn Contact Point", "Velocity","(in/s)",25.4e-3,"(deg)",np.pi/180.0)

    (xducer_displ_trange,xducer_displspec_frange,xducer_displ_timedomain,xducer_displspec,xducer_displ_filtered_timedomain,xducer_displspec_filtered) = read_spectrum_ansys(seg_table_names,seg_tables,dt,impulseexcitation_width,"Horn Contact Point", "Deflection","(in)",25.4e-3,"(deg)",np.pi/180.0)


    if laser_name is not None:
        (laser_vel_trange,laser_velspec_frange,laser_vel_timedomain,laser_velspec,laser_vel_filtered_timedomain,laser_velspec_filtered)=read_spectrum_ansys(seg_table_names,seg_tables,dt,impulseexcitation_width,laser_name, "Velocity","(in/s)",25.4e-3,"(deg)",np.pi/180.0)
        
        (laser_displ_trange,laser_displspec_frange,laser_displ_timedomain,laser_displspec,laser_displ_filtered_timedomain,laser_displspec_filtered) = read_spectrum_ansys(seg_table_names,seg_tables,dt,impulseexcitation_width,laser_name, "Deflection","(in)",25.4e-3,"(deg)",np.pi/180.0)
        
        pass

    if crack_name is not None:
        (crackcenternormalstress_trange,crackcenternormalstressspec_frange,crackcenternormalstress_timedomain,crackcenternormalstressspec,crackcenternormalstress_filtered_timedomain,crackcenternormalstressspec_filtered)=read_spectrum_ansys(seg_table_names,seg_tables,dt,impulseexcitation_width,crack_name, "Normal Stress","(psi)",6894.76,"(deg)",np.pi/180.0)
        
        
        (crackcentershearstressmajor_trange,crackcentershearstressmajorspec_frange,crackcentershearstressmajor_timedomain,crackcentershearstressmajorspec,crackcentershearstressmajor_filtered_timedomain,crackcentershearstressmajorspec_filtered)=read_spectrum_ansys(seg_table_names,seg_tables,dt,impulseexcitation_width,crack_name, "Shear Stress","(psi)",6894.76,"(deg)",np.pi/180.0)

        # !!!**** Show should shear along minor crack axis be labelled? 

        pass



    plot_paths=collections.OrderedDict()

    if plotdir is not None:
        from matplotlib import pyplot as pl                       
        pl.figure()
        # This plot should match timedomain velocity.
        # the minus sign is because the timedomain model calculation
        # depending on which is used, evaluates positive in the
        # outward normal direction.
        #if xducer_veltime_filepath is not None:
        #    xducer_veltime_args=xducer_veltime_trange,-xducer_veltime,'-'
        #    pass
        #else:
        xducer_veltime_args=[]
        #    pass
        #
        pl.plot(xducer_vel_trange,np.real(xducer_vel_filtered_timedomain),'-',
                xducer_vel_trange,np.real(xducer_vel_timedomain),'-',
                *xducer_veltime_args)
        pl.title('Transducer velocity (time domain)')
        pl.xlabel('Time(s)')
        pl.ylabel('(m/s)/(N*s)')
        pl.grid(True)
        pl.legend(('FD fusion filtered','FD fusion unfiltered','Time domain'))
        transducer_velocity_time_domain_plotpath = os.path.join(plotdir,plotprefix+'transducer_velocity_time_domain.png')
        pl.savefig(transducer_velocity_time_domain_plotpath,dpi=300)
        plot_paths["transducer_velocity_time_domain"]=transducer_velocity_time_domain_plotpath

        pl.figure()
        pl.loglog(xducer_velspec_frange,np.abs(xducer_velspec_filtered),'-',
                  xducer_velspec_frange,np.abs(xducer_velspec),'-')
        pl.axis((10,1e5,np.abs(xducer_velspec_filtered[1]),np.max(np.abs(xducer_velspec_filtered))*1.2))
        pl.xlabel('Frequency (Hz)')
        pl.ylabel('(m/s)/N')
        pl.title('Transducer velocity (frequency domain)')
        pl.grid(True)
        transducer_velocity_frequency_domain_plotpath = os.path.join(plotdir,plotprefix+'transducer_velocity_frequency_domain.png')
        pl.savefig(transducer_velocity_frequency_domain_plotpath,dpi=300)
        plot_paths["transducer_velocity_frequency_domain"]=transducer_velocity_frequency_domain_plotpath

        
        pl.figure()
        pl.semilogx(xducer_velspec_frange,np.angle(xducer_velspec_filtered))
        pl.grid(True)
        pl.xlabel('Frequency (Hz)')
        pl.ylabel('Radians')
        pl.title('Transducer phase (frequency domain)')
        transducer_phase_frequency_domain_plotpath = os.path.join(plotdir,plotprefix+'transducer_phase_frequency_domain.png')
        pl.savefig(transducer_phase_frequency_domain_plotpath,dpi=300)
        plot_paths["transducer_phase_frequency_domain"]=transducer_phase_frequency_domain_plotpath
        

        pl.figure()
        # This plot should match timedomain velocity.
        # the minus sign is because the timedomain model calculation
        # depending on which is used, evaluates positive in the
        # outward normal direction.
        #if xducer_displtime_filepath is not None:
        #    xducer_displtime_args=[xducer_displtime_trange,-xducer_displtime,'-']
        #    pass
        #else:
        xducer_displtime_args=[]
        #    pass
            
        pl.plot(xducer_displ_trange,np.real(xducer_displ_filtered_timedomain),'-',
                xducer_displ_trange,np.real(xducer_displ_timedomain),'-',
                *xducer_displtime_args)
        pl.title('Transducer displacement (time domain)')
        pl.xlabel('Time(s)')
        pl.ylabel('m/(N*s)')
        pl.grid(True)
        pl.legend(('FD fusion filtered','FD fusion unfiltered','Time domain'))
        transducer_displacement_time_domain_plotpath = os.path.join(plotdir,plotprefix+'transducer_displacement_time_domain.png')
        pl.savefig(transducer_displacement_time_domain_plotpath,dpi=300)
        plot_paths["transducer_displacement_time_domain"]=transducer_displacement_time_domain_plotpath
        
        
        if laser_name is not None:
            pl.figure()
            # This plot should match timedomain velocity.
            # the minus sign is because the timedomain model calculation
            # depending on which is used, evaluates positive in the
            # outward normal direction.
            #if laser_veltime_filepath is not None:
            #    laser_veltime_args=[laser_veltime_trange,-laser_veltime,'-']
            #    pass
            #else:
            laser_veltime_args=[]
            #    pass
            pl.plot(laser_vel_trange,np.real(laser_vel_filtered_timedomain),'-',
                    laser_vel_trange,np.real(laser_vel_timedomain),'-',
                    *laser_veltime_args)
            pl.title('Laser velocity (time domain)')
            pl.xlabel('Time(s)')
            pl.ylabel('(m/s)/(N*s)')
            pl.grid(True)
            pl.legend(('FD fusion filtered','FD fusion unfiltered','Time domain'))
            laser_velocity_time_domain_plotpath = os.path.join(plotdir,plotprefix+'laser_velocity_time_domain.png')
            pl.savefig(laser_velocity_time_domain_plotpath,dpi=300)
            plot_paths["laser_velocity_time_domain"]=laser_velocity_time_domain_plotpath


            pl.figure()
            pl.loglog(laser_velspec_frange,np.abs(laser_velspec_filtered),'-',
                      laser_velspec_frange,np.abs(laser_velspec),'-')
            pl.axis((10,1e5,np.abs(laser_velspec_filtered[1]),np.max(np.abs(laser_velspec_filtered))*1.2))
            pl.xlabel('Frequency (Hz)')
            pl.ylabel('(m/s)/N')
            pl.title('Laser velocity (frequency domain)')
            pl.grid(True)
            laser_velocity_frequency_domain_plotpath = os.path.join(plotdir,plotprefix+'laser_velocity_frequency_domain.png')
            pl.savefig(laser_velocity_frequency_domain_plotpath,dpi=300)
            plot_paths["laser_velocity_frequency_domain"]=laser_velocity_frequency_domain_plotpath
            
        
            pl.figure()
            pl.semilogx(laser_velspec_frange,np.angle(laser_velspec_filtered))
            pl.grid(True)
            pl.xlabel('Frequency (Hz)')
            pl.ylabel('Radians')
            pl.title('Laser phase (frequency domain)')
            laser_phase_frequency_domain_plotpath = os.path.join(plotdir,plotprefix+'laser_phase_frequency_domain.png')
            pl.savefig(laser_phase_frequency_domain_plotpath,dpi=300)
            plot_paths["laser_phase_frequency_domain"]=laser_phase_frequency_domain_plotpath
            pass # laser_velspec_filepaths
        
        if laser_name is not None:
            pl.figure()
            # This plot should match timedomain velocity.
            # the minus sign is because the timedomain model calculation
            # depending on which is used, evaluates positive in the
            # outward normal direction.
            #if laser_displtime_filepath is not None:
            #    laser_displtime_args=[laser_displtime_trange,-laser_displtime,'-']
            #    pass
            #else:
            laser_displtime_args=[]
            #    pass
                
            pl.plot(laser_displ_trange,np.real(laser_displ_filtered_timedomain),'-',
                    laser_displ_trange,np.real(laser_displ_timedomain),'-',
                    *laser_displtime_args)
            pl.title('Laser displacement (time domain)')
            pl.xlabel('Time(s)')
            pl.ylabel('m/(N*s)')
            pl.grid(True)
            pl.legend(('FD fusion filtered','FD fusion unfiltered','Time domain'))
            laser_displacement_time_domain_plotpath = os.path.join(plotdir,plotprefix+'laser_displacement_time_domain.png')
            pl.savefig(laser_displacement_time_domain_plotpath,dpi=300)
            plot_paths["laser_displacement_time_domain"]=laser_displacement_time_domain_plotpath
            pass
        
        
        
        if crack_name is not None:
            pl.figure()
            pl.plot(crackcenternormalstress_trange,np.real(crackcenternormalstress_filtered_timedomain),'-',
                    crackcenternormalstress_trange,np.real(crackcenternormalstress_timedomain),'-')
            pl.title('Crack center normal stress (time domain)')
            pl.xlabel('Time(s)')
            pl.ylabel('Pa/(N*s)')
            pl.grid(True)
            pl.legend(('FD fusion filtered','FD fusion unfiltered'))
            crack_center_normal_stress_time_domain_plotpath = os.path.join(plotdir,plotprefix+'crack_center_normal_stress_time_domain.png')
            pl.savefig(crack_center_normal_stress_time_domain_plotpath,dpi=300)
            plot_paths["crack_center_normal_stress_time_domain"]=crack_center_normal_stress_time_domain_plotpath
            
            
            pl.figure()
            pl.plot(crackcentershearstressmajor_trange,np.real(crackcentershearstressmajor_filtered_timedomain),'-',
                    crackcentershearstressmajor_trange,np.real(crackcentershearstressmajor_timedomain),'-')
            pl.title('Crack center shear stress (time domain)')
            pl.xlabel('Time(s)')
            pl.ylabel('Pa/(N*s)')
            pl.grid(True)
            pl.legend(('FD fusion filtered','FD fusion unfiltered'))
            crack_center_shear_stress_major_time_domain_plotpath = os.path.join(plotdir,plotprefix+'crack_center_shear_stress_major_time_domain.png')
            pl.savefig(crack_center_shear_stress_major_time_domain_plotpath,dpi=300)
            plot_paths["crack_center_shear_stress_major_time_domain"]=crack_center_shear_stress_major_time_domain_plotpath


            
            
            
            
            pl.figure()
            pl.loglog(crackcenternormalstressspec_frange,np.abs(crackcenternormalstressspec_filtered),'-',
                      crackcenternormalstressspec_frange,np.abs(crackcenternormalstressspec),'-')
            pl.axis((10,1e5,np.abs(crackcenternormalstressspec_filtered[1]),np.max(np.abs(crackcenternormalstressspec_filtered))*1.2))
            pl.xlabel('Frequency (Hz)')
            pl.ylabel('Pa/N')
            pl.title('Crack center normal stress (frequency domain)')
            pl.grid(True)
            crack_center_normal_stress_frequency_domain_plotpath = os.path.join(plotdir,plotprefix+'crack_center_normal_stress_frequency_domain.png')
            pl.savefig(crack_center_normal_stress_frequency_domain_plotpath,dpi=300)
            plot_paths["crack_center_normal_stress_frequency_domain"]=crack_center_normal_stress_frequency_domain_plotpath
            
            
            pl.figure()
            pl.loglog(crackcentershearstressmajorspec_frange,np.abs(crackcentershearstressmajorspec_filtered),'-',
                      crackcentershearstressmajorspec_frange,np.abs(crackcentershearstressmajorspec),'-')
            pl.axis((10,1e5,np.abs(crackcentershearstressmajorspec_filtered[1]),np.max(np.abs(crackcentershearstressmajorspec_filtered))*1.2))
            pl.xlabel('Frequency (Hz)')
            pl.ylabel('Pa/N')
            pl.title('Crack center shear stress (frequency domain)')
            pl.grid(True)
            crack_center_shear_stress_frequency_domain_plotpath = os.path.join(plotdir,plotprefix+'crack_center_shear_stress_frequency_domain.png')
            pl.savefig(crack_center_shear_stress_frequency_domain_plotpath,dpi=300)
            plot_paths["crack_center_shear_stress_frequency_domain"]=crack_center_shear_stress_frequency_domain_plotpath
            pass # crackcenterstressspec_filepaths

        pass
        
    print("Writing %s..." % (output_filename))

    # Write output as compressed .csv
    endpoint = xducer_displ_timedomain.shape[0]-int(endcrop/dt)
    out_frame=pd.DataFrame(index=pd.Float64Index(data=np.arange(endpoint,dtype='d')*dt,dtype='d',name="Time(s)"))
    
    # specimen_resp is displacement response in m per unit N*s impulse
    # applied at the transducer location 
    out_frame.insert(len(out_frame.columns),"specimen_resp(m/(N*s))",xducer_displ_filtered_timedomain.real[:endpoint])
    
    # specimen_mobility is velocity response in m/s per unit N*s impulse
    # applied at the transducer location 
    out_frame.insert(len(out_frame.columns),"specimen_mobility(m/(N*s^2))",xducer_vel_timedomain.real[:endpoint])
    
    if laser_name is not None:
        # specimen_laser is velocity response at laser position in 
        # m/s per unit N*s impulse
        # applied at the transducer location 
        out_frame.insert(len(out_frame.columns),"specimen_laser(m/(N*s^2))",laser_vel_timedomain.real[:endpoint])
        pass

    if crack_name is not None:

        # specimen_crackcenternormalstress is normal stress response at crack center position in 
        # Pascals per unit N*s impulse
        # applied at the transducer location 
        out_frame.insert(len(out_frame.columns),"specimen_crackcenternormalstress(Pa/(N*s))",crackcenternormalstress_timedomain.real[:endpoint])


        # specimen_crackcentershearstress is normal stress response at crack center position in 
        # Pascals per unit N*s impulse
        # applied at the transducer location 
        out_frame.insert(len(out_frame.columns),"specimen_crackcentershearstressmajor(Pa/(N*s))",crackcentershearstressmajor_timedomain.real[:endpoint])
        pass
    
    
    if output_filename.endswith(".bz2"):
        out_frame.to_csv(output_filename,compression='bz2')
        pass
    else:
        out_frame.to_csv(output_filename)
        pass


    #pl.show()
    return (output_filename,plot_paths)



def read_spectrum(segment_filepaths,colnum,dt,impulseexcitation_width,comsol_descr,descr):
    """Read COMSOL spectrum... Colnum indexed starting with 0... so colnum==0 should be frequency, colnum==1 first data column, etc."""
    
    num_segments=len(segment_filepaths)

    min_df = np.inf

    freq_key='freq (Hz)'
    measurement_key = None

    #filenames = [ filename_template % (segnum+1) for segnum in range(num_segments) ]

    probe_text_data = [ read_comsol_probe_txt.read(filename) for filename in segment_filepaths ] 

    frequency_steps = [ probe_text_data_el[2][1][freq_key]-probe_text_data_el[2][0][freq_key] for probe_text_data_el in probe_text_data ]
    min_df = np.min(frequency_steps)

                        
    FreqResp = collections.OrderedDict()
    
    for segnum in range(num_segments):
        #filename = filename_template % (segnum+1)

        #(metadata,fieldheaderdata,fieldrowdata)=read_comsol_probe_txt.read(filename)
        
        (metadata,fieldheaderdata,fieldrowdata) = probe_text_data[segnum]

        #if measurement_key is None:
        measurement_key = list(fieldheaderdata.keys())[colnum]
        if measurement_key != comsol_descr:
            raise ValueError("process_multisweep on %s: Expected column #%d marked %s; got column marked %s" % (segment_filepaths[segnum],colnum,comsol_descr,measurement_key)) 
        #    print("Got measurement field %s" % (measurement_key))
        #    pass

        this_df = fieldrowdata[1][freq_key]-fieldrowdata[0][freq_key]
        #if this_df < min_df:
        #    min_df = this_df
        #    pass

        this_segdata = []

        for row in fieldrowdata:
            FreqResp[row[freq_key]]=row[measurement_key]
            this_segdata.append(row[measurement_key])
            pass

        this_segarray=np.array(this_segdata,dtype=np.complex128)

        # plot up the separate inverse transform of this segment alone, assuming all other segments are zero
        # df = 1/(n*dt)
        # dt = 1/(n*df)
        # t_total = n*dt = 1/df
        
        # interpolate to a finer dt given_dt = (1/(n_new*df)) by padding with zeros in freq domain. 

        # this_dt_uninterpolated = 1.0/(len(this_segdata)*this_df)

        # given_dt = 1.0/(n_new*this_df)
        # ->  n_new = 1.0/(given_dt*this_df)
        n_new = int(np.round(1.0/(dt*this_df)))
        
        this_freqbase = fieldrowdata[0][freq_key] + np.arange(n_new)*this_df
        
        this_timebase = np.arange(n_new)*(1.0/(n_new*this_df))

        # multiply inverse transform by n_new (to compensate for built-in normalization of ifft) 
        # by df -- quantity integrated over -- and by exponential that compensates for the frequency 
        # corresponding to our first step
        this_segtimedata = np.fft.ifft(this_segarray,n=n_new)*n_new*this_df*np.exp(2.0*np.pi*(0+1j)*this_freqbase[0]*this_timebase)

        pl.figure()
        pl.plot(this_timebase,np.real(this_segtimedata),'-')
        pl.xlabel('Time (s)')
        pl.ylabel(descr)
        pl.title('seg%d' % (segnum+1))

        pass
    
    (trange,frange,meas_timedomain,meas,meas_filtered_timedomain,meas_filtered) = interpolate_spectrum(min_df,dt,FreqResp,impulseexcitation_width)
    


    return (trange,frange,meas_timedomain,meas,meas_filtered_timedomain,meas_filtered)

def read_timedomain(filename):
    (metadata,fieldheaderdata,fieldrowdata)=read_comsol_probe_txt.read(filename)
    time_key='Time (s)'
    measurement_key = list(fieldheaderdata.keys())[1]
    this_dt = fieldrowdata[1][time_key]-fieldrowdata[0][time_key]
    
    measurement = np.array([ row[measurement_key] for row in fieldrowdata ],dtype='d')
    trange = fieldrowdata[0][time_key] + np.arange(measurement.shape[0],dtype='d')*this_dt
    return (trange,measurement)

example_code = r"""
path='.'

xducer_velspec_filename_template="cantilever_model_xducercontactprobe_velspec_seg%d.txt"
xducer_displspec_filename_template = "cantilever_model_xducercontactprobe_displspec_seg%d.txt"

# veltime/disptime files are optional (set to None otherwise)
#xducer_veltime_filename = "cantilever_model_xducercontactprobe_vel.txt"
#xducer_displtime_filename = "cantilever_model_xducercontactprobe_displ.txt"
xducer_veltime_filename = None
xducer_displtime_filename = None

laser_velspec_filename_template="cantilever_model_laser_velspec_seg%d.txt"
laser_displspec_filename_template = "cantilever_model_laser_displspec_seg%d.txt"

# veltime/disptime files are optional (set to None otherwise)
#laser_veltime_filename = "cantilever_model_laser_vel.txt"
#laser_displtime_filename = "cantilever_model_laser_displ.txt"
laser_veltime_filename = None
laser_displtime_filename = None


plotdir = "/tmp"
output_filename=os.path.join(plotdir,"specimen_modeling_out.csv.bz2")  # should be .csv.bz2 extension
"""
def process_multisweep(path,
                       xducer_velspec_filename_template,
                       xducer_displspec_filename_template,
                       laser_velspec_filename_template,
                       laser_displspec_filename_template,
                       crackcenterstressspec_filename_template,
                       output_filename,
                       xducer_veltime_filename=None,
                       xducer_displtime_filename=None,
                       laser_veltime_filename=None,
                       laser_displtime_filename=None,
                       plotdir=None,
                       plotprefix="",
                       num_segments=4,
                       dt=1e-6,
                       impulseexcitation_width=2.0e-6, # (s) COMSOL parameter -- width of Gaussian pulse
                       endcrop=250e-3 # Amount of time (s) to crop off end of unfiltered timedomain to remove anticausal portion of signal
):  # !!!*** Remember to change process_multisweep_from_files() default params with these default params ***!!!

    xducer_velspec_filepaths=[ os.path.join(path,xducer_velspec_filename_template % (segnum+1)) for segnum in range(num_segments) ]
    xducer_displspec_filepaths=[ os.path.join(path,xducer_displspec_filename_template % (segnum+1)) for segnum in range(num_segments) ]


    if xducer_veltime_filename is not None:
        xducer_veltime_filepath = os.path.join(path,xducer_veltime_filename)
        pass
    else:
        xducer_veltime_filepath = None
        pass

    if xducer_displtime_filename is not None:
        xducer_displtime_filepath = os.path.join(path,xducer_displtime_filename)
        pass
    else:
        xducer_displtime_filepath = None
        pass

    
    if laser_velspec_filename_template is not None:
        laser_velspec_filepaths = [ os.path.join(path,laser_velspec_filename_template % (segnum+1)) for segnum in range(num_segments) ]
        pass
    else:
        laser_velspec_filepaths=None
        pass

    if laser_displspec_filename_template is not None:
        laser_displspec_filepaths = [ os.path.join(path,laser_displspec_filename_template % (segnum+1)) for segnum in range(num_segments) ]
        pass
    else:
        laser_displspec_filepaths=None
        pass

        
    if laser_veltime_filename is not None:
        laser_veltime_filepath = os.path.join(path,laser_veltime_filename)
        pass
    else:
        laser_veltime_filepath = None
        pass

    if laser_displtime_filename is not None:
        laser_displtime_filepath = os.path.join(path,laser_displtime_filename)
        pass
    else:
        laser_displtime_filepath = None
        pass

    if crackcenterstressspec_filename_template is not None:
        crackcenterstressspec_filepaths = [ os.path.join(path,crackcenterstressspec_filename_template % (segnum+1)) for segnum in range(num_segments) ]
        pass
    else:
        crackcenterstressspec_filepaths = None
        pass

    return process_multisweep_from_files(xducer_velspec_filepaths,
                                         xducer_displspec_filepaths,
                                         laser_velspec_filepaths,
                                         laser_displspec_filepaths,
                                         crackcenterstressspec_filepaths,
                                         output_filename,
                                         xducer_veltime_filepath=xducer_veltime_filepath,
                                         xducer_displtime_filepath=xducer_displtime_filepath,
                                         laser_veltime_filepath=laser_veltime_filepath,
                                         laser_displtime_filepath=laser_displtime_filepath,
                                         plotdir=plotdir,
                                         plotprefix=plotprefix,
                                         num_segments=num_segments,
                                         dt=dt,
                                         impulseexcitation_width=impulseexcitation_width, 
                                         endcrop=endcrop)

def process_multisweep_from_files(xducer_velspec_filepaths,
                                  xducer_displspec_filepaths,
                                  laser_velspec_filepaths,
                                  laser_displspec_filepaths,
                                  crackcenterstressspec_filepaths,
                                  output_filename,
                                  xducer_veltime_filepath=None,
                                  xducer_displtime_filepath=None,
                                  laser_veltime_filepath=None,
                                  laser_displtime_filepath=None,
                                  plotdir=None,
                                  plotprefix="",
                                  num_segments=4,
                                  dt=1e-6,
                                  impulseexcitation_width=2.0e-6, # (s) COMSOL parameter -- width of Gaussian pulse
                                  endcrop=250e-3 # Amount of time (s) to crop off end of unfiltered timedomain to remove anticausal portion of signal
                              ):  # !!!*** Remember to change process_multisweep() default params with these default params ***!!!
    

    #endcrop=0.0

    (xducer_vel_trange,xducer_velspec_frange,xducer_vel_timedomain,xducer_velspec,xducer_vel_filtered_timedomain,xducer_velspec_filtered)=read_spectrum(xducer_velspec_filepaths,1,dt,impulseexcitation_width,"Fillmein_in_process_multisweep.py","transducer velocity (m/s)/(N*s)")

    (xducer_displ_trange,xducer_displspec_frange,xducer_displ_timedomain,xducer_displspec,xducer_displ_filtered_timedomain,xducer_displspec_filtered)=read_spectrum(xducer_displspec_filepaths,1,dt,impulseexcitation_width,"Fillmein_in_process_multisweep.py","transducer displacement m/(N*s)")

    # Try adding step to xducer_displ_timedomain
    #xducer_displ_timedomain -= .31e-3
    #xducer_displ_timedomain=np.cumsum(xducer_vel_timedomain)*dt
    #xducer_displ_timedomain[0]=.6e-3
    

    if xducer_veltime_filepath is not None:
        (xducer_veltime_trange,xducer_veltime) = read_timedomain(xducer_veltime_filepath)
        pass

        
    if xducer_displtime_filepath is not None:
        (xducer_displtime_trange,xducer_displtime) = read_timedomain(xducer_displtime_filepath)
        pass



    if laser_velspec_filepaths is not None:
        (laser_vel_trange,laser_velspec_frange,laser_vel_timedomain,laser_velspec,laser_vel_filtered_timedomain,laser_velspec_filtered)=read_spectrum(laser_velspec_filepaths,1,dt,impulseexcitation_width,"Fillmein_in_process_multisweep.py","laser point velocity (m/s)/(N*s)")
        pass

    if laser_displspec_filepaths is not None:
        (laser_displ_trange,laser_displspec_frange,laser_displ_timedomain,laser_displspec,laser_displ_filtered_timedomain,laser_displspec_filtered)=read_spectrum(laser_displspec_filepaths,1,dt,impulseexcitation_width,"Fillmein_in_process_multisweep.py","laser point displacement m/(N*s)")
        pass

    if crackcenterstressspec_filepaths is not None:
        (crackcenternormalstress_trange,crackcenternormalstressspec_frange,crackcenternormalstress_timedomain,crackcenternormalstressspec,crackcenternormalstress_filtered_timedomain,crackcenternormalstressspec_filtered)=read_spectrum(crackcenterstressspec_filepaths,2,dt,impulseexcitation_width,"Fillmein_in_process_multisweep.py","crack center normal stress (Pa)") # colnum==1 is stress magnitude, colnum==2 is normal stress, colnum==3 is shear stress major, colnum==4 is shear stress minor
        
        (crackcentershearstressmajor_trange,crackcentershearstressmajorspec_frange,crackcentershearstressmajor_timedomain,crackcentershearstressmajorspec,crackcentershearstressmajor_filtered_timedomain,crackcentershearstressmajorspec_filtered)=read_spectrum(crackcenterstressmajorspec_filepaths,3,dt,impulseexcitation_width,"Fillmein_in_process_multisweep.py","crack center shear stress major axis (Pa)")

        (crackcentershearstressminor_trange,crackcentershearstressminorspec_frange,crackcentershearstressminor_timedomain,crackcentershearstressminorspec,crackcentershearstressminor_filtered_timedomain,crackcentershearstressminorspec_filtered)=read_spectrum(crackcenterstressminorspec_filepaths,4,dt,impulseexcitation_width,"Fillmein_in_process_multisweep.py","crack center shear stress minor axis (Pa)")
        pass

    if laser_veltime_filepath is not None:
        (laser_veltime_trange,laser_veltime) = read_timedomain(laser_veltime_filepath)
        pass
        
    if laser_displtime_filepath is not None:
        (laser_displtime_trange,laser_displtime) = read_timedomain(laser_displtime_filepath)
        pass

    plot_paths=collections.OrderedDict()

    if plotdir is not None:
        from matplotlib import pyplot as pl                       
        pl.figure()
        # This plot should match COMSOL timedomain velocity.
        # the minus sign is because the COMSOL model calculation
        # depending on which is used, evaluates positive in the
        # outward normal direction.
        if xducer_veltime_filepath is not None:
            xducer_veltime_args=xducer_veltime_trange,-xducer_veltime,'-'
            pass
        else:
            xducer_veltime_args=[]
            pass
        
        pl.plot(xducer_vel_trange,np.real(xducer_vel_filtered_timedomain),'-',
                xducer_vel_trange,np.real(xducer_vel_timedomain),'-',
                *xducer_veltime_args)
        pl.title('Transducer velocity (time domain)')
        pl.xlabel('Time(s)')
        pl.ylabel('(m/s)/(N*s)')
        pl.grid(True)
        pl.legend(('FD fusion filtered','FD fusion unfiltered','Time domain'))
        transducer_velocity_time_domain_plotpath = os.path.join(plotdir,plotprefix+'transducer_velocity_time_domain.png')
        pl.savefig(transducer_velocity_time_domain_plotpath,dpi=300)
        plot_paths["transducer_velocity_time_domain"]=transducer_velocity_time_domain_plotpath

        pl.figure()
        pl.loglog(xducer_velspec_frange,np.abs(xducer_velspec_filtered),'-',
                  xducer_velspec_frange,np.abs(xducer_velspec),'-')
        pl.axis((10,1e5,np.abs(xducer_velspec_filtered[1]),np.max(np.abs(xducer_velspec_filtered))*1.2))
        pl.xlabel('Frequency (Hz)')
        pl.ylabel('(m/s)/N')
        pl.title('Transducer velocity (frequency domain)')
        pl.grid(True)
        transducer_velocity_frequency_domain_plotpath = os.path.join(plotdir,plotprefix+'transducer_velocity_frequency_domain.png')
        pl.savefig(transducer_velocity_frequency_domain_plotpath,dpi=300)
        plot_paths["transducer_velocity_frequency_domain"]=transducer_velocity_frequency_domain_plotpath

        
        pl.figure()
        pl.semilogx(xducer_velspec_frange,np.angle(xducer_velspec_filtered))
        pl.grid(True)
        pl.xlabel('Frequency (Hz)')
        pl.ylabel('Radians')
        pl.title('Transducer phase (frequency domain)')
        transducer_phase_frequency_domain_plotpath = os.path.join(plotdir,plotprefix+'transducer_phase_frequency_domain.png')
        pl.savefig(transducer_phase_frequency_domain_plotpath,dpi=300)
        plot_paths["transducer_phase_frequency_domain"]=transducer_phase_frequency_domain_plotpath
        

        pl.figure()
        # This plot should match COMSOL timedomain velocity.
        # the minus sign is because the COMSOL model calculation
        # depending on which is used, evaluates positive in the
        # outward normal direction. 
        if xducer_displtime_filepath is not None:
            xducer_displtime_args=[xducer_displtime_trange,-xducer_displtime,'-']
            pass
        else:
            xducer_displtime_args=[]
            pass
            
        pl.plot(xducer_displ_trange,np.real(xducer_displ_filtered_timedomain),'-',
                xducer_displ_trange,np.real(xducer_displ_timedomain),'-',
                *xducer_displtime_args)
        pl.title('Transducer displacement (time domain)')
        pl.xlabel('Time(s)')
        pl.ylabel('m/(N*s)')
        pl.grid(True)
        pl.legend(('FD fusion filtered','FD fusion unfiltered','Time domain'))
        transducer_displacement_time_domain_plotpath = os.path.join(plotdir,plotprefix+'transducer_displacement_time_domain.png')
        pl.savefig(transducer_displacement_time_domain_plotpath,dpi=300)
        plot_paths["transducer_displacement_time_domain"]=transducer_displacement_time_domain_plotpath
        
        
        if laser_velspec_filepaths is not None:
            pl.figure()
            # This plot should match COMSOL timedomain velocity.
            # the minus sign is because the COMSOL model calculation
            # depending on which is used, evaluates positive in the
            # outward normal direction. 
            if laser_veltime_filepath is not None:
                laser_veltime_args=[laser_veltime_trange,-laser_veltime,'-']
                pass
            else:
                laser_veltime_args=[]
                pass
            pl.plot(laser_vel_trange,np.real(laser_vel_filtered_timedomain),'-',
                    laser_vel_trange,np.real(laser_vel_timedomain),'-',
                    *laser_veltime_args)
            pl.title('Laser velocity (time domain)')
            pl.xlabel('Time(s)')
            pl.ylabel('(m/s)/(N*s)')
            pl.grid(True)
            pl.legend(('FD fusion filtered','FD fusion unfiltered','Time domain'))
            laser_velocity_time_domain_plotpath = os.path.join(plotdir,plotprefix+'laser_velocity_time_domain.png')
            pl.savefig(laser_velocity_time_domain_plotpath,dpi=300)
            plot_paths["laser_velocity_time_domain"]=laser_velocity_time_domain_plotpath


            pl.figure()
            pl.loglog(laser_velspec_frange,np.abs(laser_velspec_filtered),'-',
                      laser_velspec_frange,np.abs(laser_velspec),'-')
            pl.axis((10,1e5,np.abs(laser_velspec_filtered[1]),np.max(np.abs(laser_velspec_filtered))*1.2))
            pl.xlabel('Frequency (Hz)')
            pl.ylabel('(m/s)/N')
            pl.title('Laser velocity (frequency domain)')
            pl.grid(True)
            laser_velocity_frequency_domain_plotpath = os.path.join(plotdir,plotprefix+'laser_velocity_frequency_domain.png')
            pl.savefig(laser_velocity_frequency_domain_plotpath,dpi=300)
            plot_paths["laser_velocity_frequency_domain"]=laser_velocity_frequency_domain_plotpath
            
        
            pl.figure()
            pl.semilogx(laser_velspec_frange,np.angle(laser_velspec_filtered))
            pl.grid(True)
            pl.xlabel('Frequency (Hz)')
            pl.ylabel('Radians')
            pl.title('Laser phase (frequency domain)')
            laser_phase_frequency_domain_plotpath = os.path.join(plotdir,plotprefix+'laser_phase_frequency_domain.png')
            pl.savefig(laser_phase_frequency_domain_plotpath,dpi=300)
            plot_paths["laser_phase_frequency_domain"]=laser_phase_frequency_domain_plotpath
            pass # laser_velspec_filepaths
        
        if laser_displspec_filepaths is not None:
            pl.figure()
            # This plot should match COMSOL timedomain velocity.
            # the minus sign is because the COMSOL model calculation
            # depending on which is used, evaluates positive in the
            # outward normal direction. 
            if laser_displtime_filepath is not None:
                laser_displtime_args=[laser_displtime_trange,-laser_displtime,'-']
                pass
            else:
                laser_displtime_args=[]
                pass
                
            pl.plot(laser_displ_trange,np.real(laser_displ_filtered_timedomain),'-',
                    laser_displ_trange,np.real(laser_displ_timedomain),'-',
                    *laser_displtime_args)
            pl.title('Laser displacement (time domain)')
            pl.xlabel('Time(s)')
            pl.ylabel('m/(N*s)')
            pl.grid(True)
            pl.legend(('FD fusion filtered','FD fusion unfiltered','Time domain'))
            laser_displacement_time_domain_plotpath = os.path.join(plotdir,plotprefix+'laser_displacement_time_domain.png')
            pl.savefig(laser_displacement_time_domain_plotpath,dpi=300)
            plot_paths["laser_displacement_time_domain"]=laser_displacement_time_domain_plotpath
            pass



        if crackcenterstressspec_filepaths is not None:
            pl.figure()
            pl.plot(crackcenternormalstress_trange,np.real(crackcenternormalstress_filtered_timedomain),'-',
                    crackcenternormalstress_trange,np.real(crackcenternormalstress_timedomain),'-')
            pl.title('Crack center normal stress (time domain)')
            pl.xlabel('Time(s)')
            pl.ylabel('Pa/(N*s)')
            pl.grid(True)
            pl.legend(('FD fusion filtered','FD fusion unfiltered'))
            crack_center_normal_stress_time_domain_plotpath = os.path.join(plotdir,plotprefix+'crack_center_normal_stress_time_domain.png')
            pl.savefig(crack_center_normal_stress_time_domain_plotpath,dpi=300)
            plot_paths["crack_center_normal_stress_time_domain"]=crack_center_normal_stress_time_domain_plotpath
            
            
            pl.figure()
            pl.plot(crackcentershearstressmajor_trange,np.real(crackcentershearstressmajor_filtered_timedomain),'-',
                    crackcentershearstressmajor_trange,np.real(crackcentershearstressmajor_timedomain),'-')
            pl.title('Crack center shear stress major (time domain)')
            pl.xlabel('Time(s)')
            pl.ylabel('Pa/(N*s)')
            pl.grid(True)
            pl.legend(('FD fusion filtered','FD fusion unfiltered'))
            crack_center_shear_stress_major_time_domain_plotpath = os.path.join(plotdir,plotprefix+'crack_center_shear_stress_major time_domain.png')
            pl.savefig(crack_center_shear_stress_major_time_domain_plotpath,dpi=300)
            plot_paths["crack_center_shear_stress_major_time_domain"]=crack_center_shear_stress_major_time_domain_plotpath

            pl.figure()
            pl.plot(crackcentershearstressminor_trange,np.real(crackcentershearstressminor_filtered_timedomain),'-',
                    crackcentershearstressminor_trange,np.real(crackcentershearstressminor_timedomain),'-')
            pl.title('Crack center shear stress minor (time domain)')
            pl.xlabel('Time(s)')
            pl.ylabel('Pa/(N*s)')
            pl.grid(True)
            pl.legend(('FD fusion filtered','FD fusion unfiltered'))
            crack_center_shear_stress_minor_time_domain_plotpath = os.path.join(plotdir,plotprefix+'crack_center_shear_stress_minor_time_domain.png')
            pl.savefig(crack_center_shear_stress_minor_time_domain_plotpath,dpi=300)
            plot_paths["crack_center_shear_stress_minor_time_domain"]=crack_center_shear_stress_minor_time_domain_plotpath
            
            
            
            
            pl.figure()
            pl.loglog(crackcenternormalstressspec_frange,np.abs(crackcenternormalstressspec_filtered),'-',
                      crackcenternormalstressspec_frange,np.abs(crackcenternormalstressspec),'-')
            pl.axis((10,1e5,np.abs(crackcenternormalstressspec_filtered[1]),np.max(np.abs(crackcenternormalstressspec_filtered))*1.2))
            pl.xlabel('Frequency (Hz)')
            pl.ylabel('Pa/N')
            pl.title('Crack center normal stress (frequency domain)')
            pl.grid(True)
            crack_center_normal_stress_frequency_domain_plotpath = os.path.join(plotdir,plotprefix+'crack_center_normal_stress_frequency_domain.png')
            pl.savefig(crack_center_normal_stress_frequency_domain_plotpath,dpi=300)
            plot_paths["crack_center_normal_stress_frequency_domain"]=crack_center_normal_stress_frequency_domain_plotpath
            
            
            pl.figure()
            pl.loglog(crackcentershearstressmajorspec_frange,np.abs(crackcentershearstressmajorspec_filtered),'-',
                      crackcentershearstressmajorspec_frange,np.abs(crackcentershearstressmajorspec),'-')
            pl.axis((10,1e5,np.abs(crackcentershearstressmajorspec_filtered[1]),np.max(np.abs(crackcentershearstressmajorspec_filtered))*1.2))
            pl.xlabel('Frequency (Hz)')
            pl.ylabel('1/N')
            pl.title('Crack center shear stress major axis (frequency domain)')
            pl.grid(True)
            crack_center_shear_stress_major_frequency_domain_plotpath = os.path.join(plotdir,plotprefix+'crack_center_shear_stress_major_frequency_domain.png')
            pl.savefig(crack_center_shear_stress_major_frequency_domain_plotpath,dpi=300)
            plot_paths["crack_center_shear_stress_major_frequency_domain"]=crack_center_shear_stress_major_frequency_domain_plotpath
            pass # crackcenterstressspec_filepaths

            
            pl.figure()
            pl.loglog(crackcentershearstressminorspec_frange,np.abs(crackcentershearstressminorspec_filtered),'-',
                      crackcentershearstressminorspec_frange,np.abs(crackcentershearstressminorspec),'-')
            pl.axis((10,1e5,np.abs(crackcentershearstressminorspec_filtered[1]),np.max(np.abs(crackcentershearstressminorspec_filtered))*1.2))
            pl.xlabel('Frequency (Hz)')
            pl.ylabel('1/N')
            pl.title('Crack center shear stress minor axis (frequency domain)')
            pl.grid(True)
            crack_center_shear_stress_minor_frequency_domain_plotpath = os.path.join(plotdir,plotprefix+'crack_center_shear_stress_minor_frequency_domain.png')
            pl.savefig(crack_center_shear_stress_minor_frequency_domain_plotpath,dpi=300)
            plot_paths["crack_center_shear_stress_minor_frequency_domain"]=crack_center_shear_stress_minor_frequency_domain_plotpath
            pass # crackcenterstressspec_filepaths

        pass
        
    print("Writing %s..." % (output_filename))

    # Write output as compressed .csv
    endpoint = xducer_displ_timedomain.shape[0]-int(endcrop/dt)
    out_frame=pd.DataFrame(index=pd.Float64Index(data=np.arange(endpoint,dtype='d')*dt,dtype='d',name="Time(s)"))
    
    # specimen_resp is displacement response in m per unit N*s impulse
    # applied at the transducer location 
    out_frame.insert(len(out_frame.columns),"specimen_resp(m/(N*s))",xducer_displ_filtered_timedomain.real[:endpoint])
    
    # specimen_mobility is velocity response in m/s per unit N*s impulse
    # applied at the transducer location 
    out_frame.insert(len(out_frame.columns),"specimen_mobility(m/(N*s^2))",xducer_vel_timedomain.real[:endpoint])
    
    if laser_velspec_filepaths is not None:
        # specimen_laser is velocity response at laser position in 
        # m/s per unit N*s impulse
        # applied at the transducer location 
        out_frame.insert(len(out_frame.columns),"specimen_laser(m/(N*s^2))",laser_vel_timedomain.real[:endpoint])
        pass

    if crackcenterstressspec_filepaths is not None:

        # specimen_crackcenternormalstress is normal stress response at crack center position in 
        # unitless per unit N*s impulse
        # applied at the transducer location 
        out_frame.insert(len(out_frame.columns),"specimen_crackcenternormalstress(Pa/(N*s))",crackcenternormalstress_timedomain.real[:endpoint])


        # specimen_crackcentershearstressmajor is shear stress response at crack center position in 
        # unitless per unit N*s impulse
        # applied at the transducer location 
        out_frame.insert(len(out_frame.columns),"specimen_crackcentershearstressmajor(Pa/(N*s))",crackcentershearstressmajor_timedomain.real[:endpoint])
        
        # specimen_crackcentershearstressminor is shear stress response at crack center position in 
        # unitless per unit N*s impulse
        # applied at the transducer location 
        out_frame.insert(len(out_frame.columns),"specimen_crackcentershearstressminor(Pa/(N*s))",crackcentershearstressminor_timedomain.real[:endpoint])
        pass
    
    
    if output_filename.endswith(".bz2"):
        out_frame.to_csv(output_filename,compression='bz2')
        pass
    else:
        out_frame.to_csv(output_filename)
        pass


    #pl.show()
    return (output_filename,plot_paths)
