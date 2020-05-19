ProcessTrak Steps
==================

.. py:function:: vibrosim_synthetic_spectrum

    This processtrak step will take the output of the modal analysis and
    calculate a spectrum. Each mode in the modal analysis is a peak in the
    spectrum, with the ambplitude and bandwidth determinted by the complex
    values eigenfrequency.

    This processtrak step is included in the ``VibroSim_Simulator`` software package.

    :param dc\:modalfreqs: |modalfreqs|
    :return dc\:modal_synthetic_spectrum: Synthetic spectrum figure.


.. py:function:: vibrocomsol_runmodal_comsol

    This processtrak step will open COMSOL, run the modal analysis study, and save the results.

    This processtrak step is included in the ``VibroSim_Simulator`` software package.

    :param dc\:dest: |dest|
    :param dc\:measident: |measident|
    :param dc\:model_comsol: |model_comsol|

    :param dc\:modalcalc_comsol: |modalcalc_comsol|
    :param dc\:modalfreqs: |modalfreqs|


.. py:function:: vibrocomsol_createdummyoutput

    This processtrak step will create a dummy heating data file. This needs to
    exist to create the COMSOL model. The processtrak step
    :py:func:`vibrocomsol_calc_heating_welder` will populate this file with
    crack heating data. The processtrak step
    :py:func:`vibrocomsol_heatflow_analysis_comsol` will use the data to analyze
    heat flows.

    This processtrak step is included in the ``VibroSim_Simulator`` software package.

    :param dc\:dest: |dest|
    :param dc\:measident: |measident|

    :return dc\:dummy_heatingdata: |dummy_heatingdata|

.. py:function:: vibrocomsol_entersweepfreqs

    This processtrak step will show the results of the modal analysis and set
    the range of frequencies to be run by a later processtrak step
    (:py:func:`vibrocomsol_sweep_analysis_comsol`). If the returns of this
    processtrak step are already defined in the experiment log (.xlp) then those
    predefined values are given as defaults for the user to accept. This would
    occur if this step has been run previously or if they are defined in the
    initial experiment log (.xlg).

    This processtrak step is included in the ``VibroSim_Simulator`` software package.

    :param dc\:modalfreqs: |modalfreqs|
    :param dc\:modalcalc_comsol: |modalcalc_comsol|

    :return dc\:sweep_start_frequency: |sweep_start_frequency|
    :return dc\:sweep_step_frequency: |sweep_step_frequency|
    :return dc\:sweep_end_frequency: |sweep_end_frequency|

.. py:function:: vibrocomsol_setsweepfrequencies_comsol

    This processtrak step will open COMSOL and set the parameters needed for
    the sweep analysis.
    
    This processtrak step is included in the ``VibroSim_Simulator`` software package.

    :param dc\:dest: |dest|
    :param dc\:measident: |measident|
    :param dc\:model_comsol: |model_comsol|
    :param dc\:sweep_start_frequency: |sweep_start_frequency|
    :param dc\:sweep_step_frequency: |sweep_step_frequency|
    :param dc\:sweep_end_frequency: |sweep_end_frequency|

    :return dc\:model_comsol_withsweepfrequencies: |model_comsol_withsweepfrequencies|

.. py:function:: vibrocomsol_sweep_analysis_comsol

    This processtrak step will open COMSOL, run the sweep analysis study, and
    save the results.

    This processtrak step is included in the ``VibroSim_Simulator`` software package.

    :param dc\:dest: |dest|
    :param dc\:measident: |measident|
    :param dc\:model_comsol_withsweepfrequencies: |model_comsol_withsweepfrequencies|

    :return dc\:sweep_spectrum: |sweep_spectrum|

.. py:function:: vibrocomsol_enterburstfreq 

    This processtrak step will show the results of the modal analysis and the
    sweep analysis, and then set the single frequency to be run by a later
    processtrak step (:py:func:`vibrocomsol_burst_analysis_comsol`). If the
    returns of this processtrak step are already defined in the experiment log
    (.xlp) then those predefined values are given as defaults for the user to
    accept. This would occur if this step has been run previously or if they
    are defined in the initial experiment log (.xlg).
    
    This processtrak step is included in the ``VibroSim_Simulator`` software package.

    :param dc\:modalfreqs: |modalfreqs|
    :param dc\:modalcalc_comsol: |modalcalc_comsol|
    :param dc\:sweep_spectrum: |sweep_spectrum|

    :return dc\:excitation_frequency: |excitation_frequency|

.. py:function:: vibrocomsol_setburstfrequency_comsol

    This processtrak step will open COMSOL and set the parameters needed for
    the burst analysis.
    
    This processtrak step is included in the ``VibroSim_Simulator`` software package.

    :param dc\:dest: |dest|
    :param dc\:measident: |measident|
    :param dc\:model_comsol: |model_comsol|
    :param dc\:excitation_frequency: |excitation_frequency|

    :return dc\:model_comsol_withburstfrequency: |model_comsol_withburstfrequency|


.. py:function:: vibrocomsol_burst_analysis_comsol

    This processtrak step will open COMSOL, run the burst analysis study, and
    save the results.

    This processtrak step is included in the ``VibroSim_Simulator`` software package.

    :param dc\:dest: |dest|
    :param dc\:measident: |measident|
    :param dc\:model_comsol_withburstfrequency: |model_comsol_withburstfrequency|

    :return dc\:burstcalc_comsol: |burstcalc_comsol|
    :return dc\:harmonicburst_normalstress: |harmonicburst_normalstress|
    :return dc\:harmonicburst_shearstressmajor: |harmonicburst_shearstressmajor|
    :return dc\:harmonicburst_shearstressminor: |harmonicburst_shearstressminor|

.. py:function:: vibrosim_calc_heating_singlefrequency

    This processtrak step will calculate the heat generated by the crack when
    the sample is subjected to the burst excitation.

    This processtrak step is included in the ``VibroSim_Simulator`` software package.

    :param dc\:dest: |dest|
    :param dc\:measident: |measident|
    :param dc\:friction_coefficient: |friction_coefficient|
    :param dc\:msqrtR: |msqrtR|
    :param dc\:staticload: |staticload|
    :param dc\:exclength: |exclength|
    :param dc\:tortuosity: |tortuosity|
    :param dc\:numdraws: |numdraws|
    :param dc\:YoungsModulus: |YoungsModulus|
    :param dc\:PoissonsRatio: |PoissonsRatio|
    :param dc\:YieldStrength: |YieldStrength|
    :param dc\:Density: |Density|
    :param dc\:crack_model_normal: |crack_model_normal|
    :param dc\:crack_model_shear: |crack_model_shear|
    :param dc\:crack_model_shear_factor: |crack_model_shear_factor|
    :param dc\:exc_t0: |exc_t0|
    :param dc\:exc_t1: |exc_t1|
    :param dc\:exc_t2: |exc_t2|
    :param dc\:exc_t3: |exc_t3|
    :param dc\:exc_t4: |exc_t4|
    :param dc\:excitation_frequency: |excitation_frequency|
    :param dc\:harmonicburst_normalstress: |harmonicburst_normalstress|
    :param dc\:harmonicburst_shearstressmajor: |harmonicburst_shearstressmajor|
    :param dc\:harmonicburst_shearstressminor: |harmonicburst_shearstressminor|
    :param dc\:crack_type_side1: |crack_type_side1|
    :param dc\:crack_type_side2: |crack_type_side2|
    :param dc\:crack_shearstress_axis: |crack_shearstress_axis|
    :param dc\:thickness: |thickness|
    :param dc\:closurestate_side1: |closurestate_side1|
    :param dc\:closurestate_side2: |closurestate_side2|
    :param dc\:a_side1: |a_side1|
    :param dc\:a_side2: |a_side2|

    :return dc\:heatpower: |heatpower|
    :return dc\:heatingdata: |heatingdata|
    :return dc\:heatingtotalpower: |heatingtotalpower|

.. py:function:: vibrocomsol_heatflow_analysis_comsol

    This processtrak step will open COMSOL, run the heatflow analysis study, and save the results.

    This processtrak step is included in the ``VibroSim_Simulator`` software package.

    :param dc\:dest: |dest|
    :param dc\:measident: |measident|
    :param dc\:model_comsol: |model_comsol|
    :param dc\:heatingdata: |heatingdata|
    :param dc\:exc_t3: |exc_t3|

    :return dc\:vibro_heating_image: |vibro_heating_image|
    :return dc\:heatflow_comsol: |heatflow_comsol|


.. py:function:: vibrosim_calc_heating_welder

    This processtrak step will calculate the heat generated by the crack when
    the sample is subjected to the welder excitation.

    :param dc\:dest: |dest|
    :param dc\:measident: |measident|
    :param dc\:friction_coefficient: |friction_coefficient|
    :param dc\:msqrtR: |msqrtR|
    :param dc\:staticload: |staticload|
    :param dc\:tortuosity: |tortuosity|
    :param dc\:numdraws: |numdraws|
    :param dc\:YoungsModulus: |YoungsModulus|
    :param dc\:PoissonsRatio: |PoissonsRatio|
    :param dc\:YieldStrength: |YieldStrength|
    :param dc\:Density: |Density|
    :param dc\:crack_model_normal: |crack_model_normal|
    :param dc\:crack_model_shear: |crack_model_shear|
    :param dc\:crack_model_shear_factor: |crack_model_shear_factor|
    :param dc\:exc_t0: |exc_t0|
    :param dc\:exc_t4: |exc_t4|
    :param dc\:motion: |motion|
    :param dc\:crack_type_side1: |crack_type_side1|
    :param dc\:crack_type_side2: |crack_type_side2|
    :param dc\:crack_shearstress_axis: |crack_shearstress_axis|
    :param dc\:thickness: |thickness|
    :param dc\:closurestate_side1: |closurestate_side1|
    :param dc\:closurestate_side2: |closurestate_side2|
    :param dc\:a_side1: |a_side1|
    :param dc\:a_side2: |a_side2|

    :return dc\:heatingdata: |heatingdata|
    :return dc\:heatingtotalpower: |heatingtotalpower|

.. py:function:: vibrosim_eval_closure_state_from_tip_positions
    
    This processtrak step is used in vibrosim to evaluate crack closure state
    from crack tip positions given in an XML element. It is provided by the
    ``crackclosuresim2`` package.
    
    The crack closure state is given as four arrays interpreted as text
    within XML elements of the experiment log,
    e.g: ::

        <dc:reff_side1 dcv:units="m" dcv:arraystorageorder="C">
          <dcv:arrayshape>9</dcv:arrayshape>
          <dcv:arraydata>
            .5e-3 .7e-3 .9e-3 1.05e-3 1.2e-3 1.33e-3 1.45e-3 1.56e-3 1.66e-3
          </dcv:arraydata>
        </dc:reff_side1>
        <dc:seff_side1 dcv:units="Pa" dcv:arraystorageorder="C">
          <dcv:arrayshape>9</dcv:arrayshape>
          <dcv:arraydata>
            0.0 50e6 100e6 150e6 200e6 250e6 300e6 350e6 400e6
          </dcv:arraydata>
        </dc:seff_side1>
        
        <dc:reff_side2 dcv:units="m" dcv:arraystorageorder="C">
          <dcv:arrayshape>9</dcv:arrayshape>
          <dcv:arraydata>
            .5e-3 .7e-3 .9e-3 1.05e-3 1.2e-3 1.33e-3 1.45e-3 1.56e-3 1.66e-3
          </dcv:arraydata>
        </dc:reff_side2>
        <dc:seff_side2 dcv:units="Pa" dcv:arraystorageorder="C">
          <dcv:arrayshape>9</dcv:arrayshape>
          <dcv:arraydata>
            0.0 50e6 100e6 150e6 200e6 250e6 300e6 350e6 400e6
          </dcv:arraydata>
        </dc:seff_side2>
    
    The ``reff`` (effective tip radius) values are given in meters and the
    ``seff`` (corresponding normal stress) values are given in Pascals.
    The radius values should be listed in increasing order. The last
    radius value on each side (side1 - left or side2 - right) should
    correspond to the length of that side of the crack. 

    :param dc\:dest: |dest|
    :param dc\:measident: |measident|
    :param dc\:YoungsModulus: |YoungsModulus|
    :param dc\:PoissonsRatio: |PoissonsRatio|
    :param dc\:YieldStrength: |YieldStrength|
    :param dc\:reff_side1: |reff_side1|
    :param dc\:seff_side1: |seff_side1|
    :param dc\:reff_side2: |reff_side2|
    :param dc\:seff_side2: |seff_side2|
    :param dc\:crack_model_normal: |crack_model_normal|
    :param dc\:crack_model_shear: |crack_model_shear|

    :return dc\:closureplot_side1: |closureplot_side1| 
    :return dc\:closureplot_side2: |closureplot_side2| 
    :return dc\:closurestate_side1: |closurestate_side1| 
    :return dc\:closurestate_side2: |closurestate_side2| 
    :return dc\:a_side1: |a_side1| 
    :return dc\:a_side2: |a_side2| 
  

.. py:function:: vibrosim_plot_welder_motion
    
    This processtrak step will plot the welder motion.

    Provided by the ``VibroSim_WelderModel`` package.
    
    :param dc\:dest: |dest|
    :param dc\:measident: |measident|
    :param dc\:motion: |motion|
    :param dc\:exc_t0: |exc_t0|

    :return plots: Welder motion plots.

.. py:function:: vibrosim_simulate_welder
    
    Provided by the ``VibroSim_WelderModel`` package.

    :param dc\:dest: |dest|
    :param dc\:measident: |measident|
    :param dc\:dynamicmodel: |dynamicmodel|
    :param dc\:exc_t0: |exc_t0|
    :param dc\:exc_t1: |exc_t1|
    :param dc\:exc_t2: |exc_t2|
    :param dc\:exc_t3: |exc_t3|
    :param dc\:exc_t4: |exc_t4|
    :param dc\:mass_of_welder_and_slider: |mass_of_welder_and_slider|
    :param dc\:pneumatic_force: |pneumatic_force|
    :param dc\:welder_elec_ampl: |welder_elec_ampl|
    :param dc\:YoungsModulus: |YoungsModulus|
    :param dc\:PoissonsRatio: |PoissonsRatio|
    :param dc\:welder_spring_constant: |welder_spring_constant|
    :param dc\:R_contact: |R_contact|
    :param dc\:welder_elec_freq: |welder_elec_freq|
    :param dc\:contact_model_timestep: |contact_model_timestep|
    :param dc\:gpu_device_priority_list: |gpu_device_priority_list|
    :param dc\:gpu_precision: |gpu_precision|

    :return dc\:motion: |motion|


.. py:function:: vibrosim_process_multisweep

    Description

    :param dc\:dest: |dest|
    :param dc\:measident: |measident|
    :param dc\:seg1_xducercontactprobe_displ: |segX_xducercontactprobe_displ|
    :param dc\:seg1_xducercontactprobe_vel: |segX_xducercontactprobe_vel|
    :param dc\:seg1_laser_displ: |segX_laser_displ|
    :param dc\:seg1_laser_vel: |segX_laser_vel|
    :param dc\:seg1_crackcenterstress: |segX_crackcenterstress|
    :param dc\:seg2_xducercontactprobe_displ: |segX_xducercontactprobe_displ|
    :param dc\:seg2_xducercontactprobe_vel: |segX_xducercontactprobe_vel|
    :param dc\:seg2_laser_displ: |segX_laser_displ|
    :param dc\:seg2_laser_vel: |segX_laser_vel|
    :param dc\:seg2_crackcenterstress: |segX_crackcenterstress|
    :param dc\:seg3_xducercontactprobe_displ: |segX_xducercontactprobe_displ|
    :param dc\:seg3_xducercontactprobe_vel: |segX_xducercontactprobe_vel|
    :param dc\:seg3_laser_displ: |segX_laser_displ|
    :param dc\:seg3_laser_vel: |segX_laser_vel|
    :param dc\:seg3_crackcenterstress: |segX_crackcenterstress|
    :param dc\:seg4_xducercontactprobe_displ: |segX_xducercontactprobe_displ|
    :param dc\:seg4_xducercontactprobe_vel: |segX_xducercontactprobe_vel|
    :param dc\:seg4_laser_displ: |segX_laser_displ|
    :param dc\:seg4_laser_vel: |segX_laser_vel|
    :param dc\:seg4_crackcenterstress: |segX_crackcenterstress|
    :param dc\:endcrop: |endcrop|
                                                 
    :return dc\:dynamicmodel: not_implemented

.. py:function:: vibrocomsol_multisweep_seg_analysis_comsol

    Description

    :param dc\:dest: |dest|
    :param dc\:measident: |measident|
    :param dc\:model_comsol_withsegboundaries: |model_comsol_withsegboundaries|
    :param segnum_int: |segnum_int|

    :return dc\:segX_xducercontactprobe_displ: |segX_xducercontactprobe_displ|
    :return dc\:segX_xducercontactprobe_vel: |segX_xducercontactprobe_vel|
    :return dc\:segX_laser_displ: |segX_laser_displ|
    :return dc\:segX_laser_vel: |segX_laser_vel|
    :return dc\:segX_crackcenterstress: |segX_crackcenterstress|

.. py:function:: vibrocomsol_set_freqbands_comsol

    This processtrak step opens a COMSOL file and sets the frequency bands for the multisweep study.

    :param dc\:dest: |dest|
    :param dc\:measident: |measident|
    :param dc\:model_comsol: |model_comsol|
    :param dc\:freqband_seg1_start: |freqband_segX_start|
    :param dc\:freqband_seg1_step: |freqband_segX_step|
    :param dc\:freqband_seg1_end: |freqband_segX_end|
    :param dc\:freqband_seg2_start: |freqband_segX_start|
    :param dc\:freqband_seg2_step: |freqband_segX_step|
    :param dc\:freqband_seg2_end: |freqband_segX_end|
    :param dc\:freqband_seg3_start: |freqband_segX_start|
    :param dc\:freqband_seg3_step: |freqband_segX_step|
    :param dc\:freqband_seg3_end: |freqband_segX_end|
    :param dc\:freqband_seg4_start: |freqband_segX_start|
    :param dc\:freqband_seg4_step: |freqband_segX_step|
    :param dc\:freqband_seg4_end: |freqband_segX_end|

    :return dc\:model_comsol_withsegboundaries: |model_comsol_withsegboundaries|

.. py:function:: vibrosim_optimize_freqbands

    This processtrak step optimizes the frequency bands for the processtrak
    step :py:func:`vibrocomsol_multisweep_seg_analysis_comsol`. 

    Run this on output of modal analysis to interpret 
    the modal decay coefficients and plan a three or four segment
    frequency domain calculation that will be invertable 
    to a time-domain response. 

    :param dc\:modalfreqs: |modalfreqs|
    :param dc\:temporal_decay_divisor: |temporal_decay_divisor|
    :param dc\:spectral_decay_divisor: |spectral_decay_divisor|

    :return  dc\:freqband_seg1_start: |freqband_segX_start|
    :return  dc\:freqband_seg1_step: |freqband_segX_step|
    :return  dc\:freqband_seg1_end: |freqband_segX_end|
    :return  dc\:freqband_seg2_start: |freqband_segX_start|
    :return  dc\:freqband_seg2_step: |freqband_segX_step|
    :return  dc\:freqband_seg2_end: |freqband_segX_end|
    :return  dc\:freqband_seg3_start: |freqband_segX_start|
    :return  dc\:freqband_seg3_step: |freqband_segX_step|
    :return  dc\:freqband_seg3_end: |freqband_segX_end|
    :return  dc\:freqband_seg4_start: |freqband_segX_start|
    :return  dc\:freqband_seg4_step: |freqband_segX_step|
    :return  dc\:freqband_seg4_end: |freqband_segX_end|


.. |modalfreqs| replace:: Frequencies of the eigenmodes.
.. |modalcalc_comsol| replace:: Save file for the COMSOL model with modal results.
.. |dest| replace:: Designated results output folder.
.. |measident| replace:: Measurement identifier.
.. |model_comsol| replace:: Save file for the COMSOL model.
.. |dummy_heatingdata| replace:: Output file for heating data.

.. |sweep_start_frequency| replace:: Starting frequency for a sweep analysis.
.. |sweep_step_frequency| replace:: Frequency step for a sweep analysis.
.. |sweep_end_frequency| replace:: Ending frequency for a sweep analysis.
.. |model_comsol_withsweepfrequencies| replace:: Save file for the COMSOL model with sweep study results.

.. |sweep_spectrum| replace:: Sweep spectrum image.

.. |excitation_frequency| replace:: Burst frequency.
.. |model_comsol_withburstfrequency| replace:: Save file for the COMSOL model with burst study results.

.. |burstcalc_comsol| replace:: Save file for the COMSOL model with burst study results.
.. |harmonicburst_normalstress| replace:: Stress in the model at the crack center, normal to the crack face. (Mode I)
.. |harmonicburst_shearstressmajor| replace:: Stress in the model at the crack center, in the crack semi-major direction. (Mode II)
.. |harmonicburst_shearstressminor| replace:: Stress in the model at the crack center, in the crack semi-minor direction. (Mode III)

.. |friction_coefficient| replace:: Friction coefficient of the crack surface.
.. |msqrtR| replace:: Crack asperity density.
.. |staticload| replace:: Static bending opening load on crack.
.. |exclength| replace:: Excitation length.
.. |tortuosity| replace:: Crack tortuosity, standard deviation of the crack trajectory.
.. |numdraws| replace:: not_implemented
.. |YoungsModulus| replace:: Youngs modulus of the material.
.. |PoissonsRatio| replace:: Poissons Ratio of the material. 
.. |YieldStrength| replace:: Yield strength of the material. 
.. |Density| replace:: Density of the material.
.. |crack_model_normal| replace:: Model used for calculating crack displacement normal to crack surface.
.. |crack_model_shear| replace:: Model used for calculating crack dispalcement in shear. 
.. |crack_model_shear_factor| replace:: Sensitivity factor for shear vs normal heating.
.. |exc_t0| replace:: Start of excitation envelope ramp-up.
.. |exc_t1| replace:: End of excitation envelope ramp-up.
.. |exc_t2| replace:: Start of excitation envelope ramp-down.
.. |exc_t3| replace:: End of excitation envelope ramp down.
.. |exc_t4| replace:: End of excitation vibration calculation.
.. |crack_type_side1| replace:: Crack type of side1, can be ``halfthrough`` or ``quarterpenny``. 
.. |crack_type_side2| replace:: Crack type of side2, can be ``halfthrough`` or ``quarterpenny``. 
.. |crack_shearstress_axis| replace:: ``major`` (mode II) or ``minor`` (mode III) axis, defines the axis used to calculate shear heating.
.. |thickness| replace:: Thickness of the material at the crack, used only for ``halfthrough`` cracks. 
.. |closurestate_side1| replace:: Closure state, side 1. 
.. |closurestate_side2| replace:: Closure state, side 2. 
.. |a_side1| replace:: Semimajor axis length of side 1.
.. |a_side2| replace:: Semimajor axis length of side 2.

.. |heatpower| replace:: Heat power vs crack location figure.
.. |heatingdata| replace:: Heating data, replaces the dummy heating data.
.. |heatingtotalpower| replace:: Total heating power of the crack.

.. |motion| replace:: not_implemented
.. |vibro_heating_image| replace:: Image of temperature data as if from IR camera.
.. |heatflow_comsol| replace:: Save file for the COMSOL model with heatflow study results.

.. |reff_side1| replace:: Effective tip radius array of crack side 1.
.. |seff_side1| replace:: Normal stress corresponding to tip radius array, side 1.
.. |reff_side2| replace:: Effective tip radius of crack side 2.
.. |seff_side2| replace:: Normal stress corresponding to tip radius array, side 2.
.. |closureplot_side1| replace:: Plot of the closure state, side 1.
.. |closureplot_side2| replace:: Plot of the closure state, side 2.

.. |dynamicmodel| replace:: not_implemented
.. |mass_of_welder_and_slider| replace:: not_implemented
.. |pneumatic_force| replace:: not_implemented
.. |welder_elec_ampl| replace:: not_implemented
.. |welder_spring_constant| replace:: not_implemented
.. |R_contact| replace:: not_implemented
.. |welder_elec_freq| replace:: not_implemented
.. |contact_model_timestep| replace:: not_implemented
.. |gpu_device_priority_list| replace:: not_implemented
.. |gpu_precision| replace:: not_implemented

.. |segX_xducercontactprobe_displ| replace:: Transducer contact probe displacement.
.. |segX_xducercontactprobe_vel| replace:: Transducer contact probe velocity.
.. |segX_laser_displ| replace:: Displacement at laser vibrometer spot. 
.. |segX_laser_vel| replace:: Velocity at laser vibrometer spot.
.. |segX_crackcenterstress| replace:: Crack center stress.
.. |endcrop| replace:: not_implemented

.. |segnum_int| replace:: not_implemented

.. |freqband_segX_start| replace:: Starting frequency of a frequency band.
.. |freqband_segX_step| replace:: Step frequency of a frequency band.
.. |freqband_segX_end| replace:: End frequency of a frequency band.
.. |model_comsol_withsegboundaries| replace:: Save file for the COMSOL model with segment boundaries.

.. |temporal_decay_divisor| replace:: not_implemented
.. |spectral_decay_divisor| replace:: not_implemented
