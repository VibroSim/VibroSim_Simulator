import sys
import os
import os.path
import numpy as np
from matplotlib import pyplot as pl

def assess_segment(dt,nt,temporal_decay_divisor,spectral_decay_divisor,complex_freqs,segstart,segentries):
    # Our requirement is that all resonances must decay to 1/1000 
    # of their peak before the FFT periodicity kicks in.
    # That is require exp(-2*pi*fi*t)=1e-3 prior to periodicity
    # or -2*pi*fi*t=ln(1e-3)
    #or t = -ln(1e-3)/(2*pi*fi) (minimum time prior to periodicity)


    df=1/(nt*dt)
    
    segfreq_idxs = (complex_freqs.real >= (segstart*df)) & (complex_freqs.imag < ((segstart+segentries)*df))
    if np.count_nonzero(segfreq_idxs) > 0:
        seg_tmaxreqd = np.max(-np.log(1.0/temporal_decay_divisor)/(2.0*np.pi*complex_freqs[segfreq_idxs].imag))
        pass
    else:
        seg_tmaxreqd = 0.0
        pass

    # The spectral response close to the resonance is 
    # 1/(i*2*pi*f-i*2*pi*fr+2*pi*fi)
    # = 1/(i*2*pi*(f-fr)+2*pi*fi)
    # When doing the next frequency band up, we don't want any resonance 
    # just below the lower cutoff to alias. So we want to consider all 
    # such modes where (fcutoff-fr) <=  10*fi in determining delta_t for
    #this new frequency band... i.e. max_(which mode) -ln(1e-3)/(2*pi*fi)
    # for all modes above the cutoff + all modes below where fcutoff-fr) <=  10*fi
    # (10 above is now spectral_decay_divisor)

    belowfreq_idxs = (complex_freqs.real < (segstart*df)) & ((segstart*df - complex_freqs.real) <= spectral_decay_divisor*complex_freqs.imag)
    if np.count_nonzero(belowfreq_idxs) > 0:
        belowfreq_tmaxreqd = np.max(-np.log(1.0/temporal_decay_divisor)/(2.0*np.pi*complex_freqs[belowfreq_idxs].imag))
        pass
    else: 
        belowfreq_tmaxreqd = 0.0
        pass

    abovefreq_idxs = (complex_freqs.real > ((segstart+segentries)*df)) & ((complex_freqs.real-(segstart+segentries)*df) <= spectral_decay_divisor*complex_freqs.imag)
    if np.count_nonzero(abovefreq_idxs):
        abovefreq_tmaxreqd = np.max(-np.log(1.0/temporal_decay_divisor)/(2.0*np.pi*complex_freqs[abovefreq_idxs].imag))
        pass
    else:
        abovefreq_tmaxreqd = 0.0
        pass
    
    net_tmaxreqd=max(seg_tmaxreqd,belowfreq_tmaxreqd,abovefreq_tmaxreqd)
    if net_tmaxreqd > 0.0:
        nt_seg = 1 << int(np.ceil(np.log(net_tmaxreqd/dt)/np.log(2)))
        segment_spectral_downsampling = nt//nt_seg
        pass
    else:
        segment_spectral_downsampling = 1024
        pass
    if segment_spectral_downsampling > 1024:
        segment_spectral_downsampling = 1024  # don't downsample more than 1024
        pass

    return segment_spectral_downsampling

def split_segment(dt,nt,temporal_decay_divisor,spectral_decay_divisor,complex_freqs,start,end):

    initial_step = int(np.ceil((end-start)/5000.0))

    cutoffrange = np.arange(start+initial_step,end,initial_step)
    totalsteps = np.zeros(cutoffrange.shape,dtype='i')


    for cutoffcnt in range(cutoffrange.shape[0]):
        cutoffentries=cutoffrange[cutoffcnt]
        
        seg1_downsampling = assess_segment(dt,nt,temporal_decay_divisor,spectral_decay_divisor,complex_freqs,start,cutoffentries)
        seg1_steps = (cutoffentries-start)//seg1_downsampling
        seg2_downsampling = assess_segment(dt,nt,temporal_decay_divisor,spectral_decay_divisor,complex_freqs,cutoffentries,end-cutoffentries)
        seg2_steps = (end-cutoffentries)//seg2_downsampling
        totalsteps[cutoffcnt] = seg1_steps+seg2_steps
        pass


    

    #pl.figure()
    #pl.plot(cutoffrange,totalsteps,'-o')
    
    bestcutoff = cutoffrange[np.argmin(totalsteps)]

    # refine bestcutoff 
    refined_cutoffrange = np.arange(max(start+1,bestcutoff-initial_step),min(end,bestcutoff+initial_step))
    refined_totalsteps = np.zeros(refined_cutoffrange.shape,dtype='i')

    for cutoffcnt in range(refined_cutoffrange.shape[0]):
        cutoffentries=refined_cutoffrange[cutoffcnt]
        
        seg1_downsampling = assess_segment(dt,nt,temporal_decay_divisor,spectral_decay_divisor,complex_freqs,start,cutoffentries)
        seg1_steps = (cutoffentries-start)//seg1_downsampling
        seg2_downsampling = assess_segment(dt,nt,temporal_decay_divisor,spectral_decay_divisor,complex_freqs,cutoffentries,end-cutoffentries)
        seg2_steps = (end-cutoffentries)//seg2_downsampling
        refined_totalsteps[cutoffcnt] = seg1_steps+seg2_steps
        pass



    seg1_downsampling = assess_segment(dt,nt,temporal_decay_divisor,spectral_decay_divisor,complex_freqs,start,bestcutoff)
    seg1_steps = (bestcutoff-start)//seg1_downsampling

    seg2_downsampling = assess_segment(dt,nt,temporal_decay_divisor,spectral_decay_divisor,complex_freqs,bestcutoff,end-bestcutoff)
    seg2_steps = (end-bestcutoff)//seg2_downsampling

    refined_bestcutoff = refined_cutoffrange[np.argmin(refined_totalsteps)]

    
    return (refined_bestcutoff,seg1_downsampling,seg1_steps,seg2_downsampling,seg2_steps)


def optimize_freqbands(filename,temporal_decay_divisor,spectral_decay_divisor,min_freq,dt,FourSegs):
    """if FourSegs is True we generate four segments, otherwise
    we generate three segments"""

    if filename.lower().endswith(".txt"):
        # .txt extension... assume COMSOL
        from VibroSim_Simulator import read_comsol_probe_txt
        (metadata,fieldheaderdata,fieldrowdata) = read_comsol_probe_txt.read(filename)
        complex_freqs_raw=fieldrowdata[0]['freq (1/s)'][1]
        pass
    elif filename.lower().endswith(".xls") or filename.lower().endswith(".xlsx"):
        # Excel... assume ANSYS
        import pandas as pd

        #pd.read_excel(filename)
        raise NotImplementedError()
        pass
    else:
        raise ValueError("Unknown filename extension %s" % (sys.path.splitext(filename)[1]))
    
    complex_freqs = complex_freqs_raw[complex_freqs_raw.real >= min_freq]
    t_maxreqd = np.max(-np.log(1.0/temporal_decay_divisor)/(2.0*np.pi*complex_freqs.imag))

    # For the moment let's assume a power of 2 number of timesteps
    nt_maxreqd = int(t_maxreqd//dt)
    nt = 1 << int(np.ceil(np.log(nt_maxreqd)/np.log(2)))
    
    print("nt_maxreqd=%d; nt=%d" % (nt_maxreqd,nt))

    if FourSegs: 
        # Split entire space
        (secondcutoff,seg1_downsampling,seg1_steps,seg2_downsampling,seg2_steps)=split_segment(dt,nt,temporal_decay_divisor,spectral_decay_divisor,complex_freqs,0,nt)

        # Split 2nd half into (seg3,seg4)
        (thirdcutoff,seg3_downsampling,seg3_steps,seg4_downsampling,seg4_steps)=split_segment(dt,nt,temporal_decay_divisor,spectral_decay_divisor,complex_freqs,secondcutoff,nt)
        
        # Split 1st half into (seg1,seg2)
        (firstcutoff,seg1_downsampling,seg1_steps,seg2_downsampling,seg2_steps)=split_segment(dt,nt,temporal_decay_divisor,spectral_decay_divisor,complex_freqs,0,secondcutoff)
        
        changes = True
        
        while changes:
            old1st = firstcutoff
            old2nd = secondcutoff
            old3rd = thirdcutoff
            
            # reoptimize second cutoff
            (secondcutoff,seg2_downsampling,seg2_steps,seg3_downsampling,seg3_steps)=split_segment(dt,nt,temporal_decay_divisor,spectral_decay_divisor,complex_freqs,firstcutoff,thirdcutoff)
            
            # reoptimize first cutoff
            (firstcutoff,seg1_downsampling,seg1_steps,seg2_downsampling,seg2_steps)=split_segment(dt,nt,temporal_decay_divisor,spectral_decay_divisor,complex_freqs,0,secondcutoff)
            
            # reoptimize third cutoff
            (thirdcutoff,seg3_downsampling,seg3_steps,seg4_downsampling,seg4_steps)=split_segment(dt,nt,temporal_decay_divisor,spectral_decay_divisor,complex_freqs,secondcutoff,nt)
            
            if old1st != firstcutoff or old2nd != secondcutoff or old3rd != thirdcutoff: 
                changes=True
                pass
            else: 
                changes=False
                pass
            pass
            
            
        print("seg1_steps=%d; seg2_steps=%d; seg3_steps=%d seg4_steps=%d" % (seg1_steps,seg2_steps,seg3_steps,seg4_steps))

        seg1_freqstart = 0
        seg1_freqstep = seg1_downsampling/(nt*dt)
        seg1_freqend = seg1_steps*seg1_downsampling/(nt*dt)

        seg2_freqstart = firstcutoff/(nt*dt)
        seg2_freqstep = seg2_downsampling/(nt*dt)
        seg2_freqend = (firstcutoff + seg2_steps*seg2_downsampling)/(nt*dt)

        seg3_freqstart = secondcutoff/(nt*dt)
        seg3_freqstep = seg3_downsampling/(nt*dt)
        seg3_freqend = (secondcutoff + seg3_steps*seg3_downsampling)/(nt*dt)

        seg4_freqstart = thirdcutoff/(nt*dt)
        seg4_freqstep = seg4_downsampling/(nt*dt)
        seg4_freqend = (nt-1)/(nt*dt)

        print("seg1 freqstart=0, freqstep = %.18g, freqend=%.18g" % (seg1_freqstep,seg1_freqend))
        print("seg2 freqstart=%.18g, freqstep = %.18g, freqend=%.18g" % (seg2_freqstart,seg2_freqstep,seg2_freqend))
        print("seg3 freqstart=%.18g, freqstep = %.18g, freqend=%.18g" % (seg3_freqstart,seg3_freqstep,seg3_freqend))
        print("seg4 freqstart=%.18g, freqstep = %.18g, freqend=%.18g" % (seg4_freqstart,seg4_freqstep,seg4_freqend))
        
        print("total_steps=%d" % (seg1_steps+seg2_steps+seg3_steps+seg4_steps))
        
        pass
    else: 
        # ThreeSegs
        
        # Split entire space
        (firstcutoff,seg1_downsampling,seg1_steps,seg2_downsampling,seg2_steps)=split_segment(dt,nt,temporal_decay_divisor,spectral_decay_divisor,complex_freqs,0,nt)
        
        # Split 2nd half into (seg2,seg3)
        (secondcutoff,seg2_downsampling,seg2_steps,seg3_downsampling,seg3_steps)=split_segment(dt,nt,temporal_decay_divisor,spectral_decay_divisor,complex_freqs,firstcutoff,nt)
        
        # reoptimize first cutoff
        (firstcutoff,seg1_downsampling,seg1_steps,seg2_downsampling,seg2_steps)=split_segment(dt,nt,temporal_decay_divisor,spectral_decay_divisor,complex_freqs,0,secondcutoff)
        
        # reoptimize second cutoff
        (secondcutoff,seg2_downsampling,seg2_steps,seg3_downsampling,seg3_steps)=split_segment(dt,nt,temporal_decay_divisor,spectral_decay_divisor,complex_freqs,firstcutoff,nt)
        
        print("seg1_steps=%d; seg2_steps=%d; seg3_steps=%d" % (seg1_steps,seg2_steps,seg3_steps))

        seg1_freqstart = 0
        seg1_freqstep = seg1_downsampling/(nt*dt)
        seg1_freqend = seg1_steps*seg1_downsampling/(nt*dt)

        seg2_freqstart = firstcutoff/(nt*dt)
        seg2_freqstep = seg2_downsampling/(nt*dt)
        seg2_freqend = (firstcutoff + seg2_steps*seg2_downsampling)/(nt*dt)

        seg3_freqstart = secondcutoff/(nt*dt)
        seg3_freqstep = seg3_downsampling/(nt*dt)
        seg3_freqend = (nt-1)/(nt*dt)

        seg4_freqstart = seg3_freqend
        seg4_freqstep = seg3_freqstep
        seg4_freqend = seg4_freqstart-seg4_freqstep

        print("seg1 freqstart=0, freqstep = %.18g, freqend=%.18g" % (seg1_freqstep,seg1_freqend))
        print("seg2 freqstart=%.18g, freqstep = %.18g, freqend=%.18g" % (seg2_freqstart,seg2_freqstep,seg2_freqend))
        print("seg3 freqstart=%.18g, freqstep = %.18g, freqend=%.18g" % (seg3_freqstart,seg3_freqstep,seg3_freqend))

        print("total_steps=%d" % (seg1_steps+seg2_steps+seg3_steps))
        pass
        
        
    return (seg1_freqstart,seg1_freqstep,seg1_freqend,
            seg2_freqstart,seg2_freqstep,seg2_freqend,
            seg3_freqstart,seg3_freqstep,seg3_freqend,
            seg4_freqstart,seg4_freqstep,seg4_freqend)
