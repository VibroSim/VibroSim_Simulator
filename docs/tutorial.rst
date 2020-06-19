Tutorial
==================

Lets look in depth at ``vibrosim_demo3`` in the examples folder. This example
has three files associated with it, ``vibrosim_demo3.prx``,
``vibrosim_demo3.xlg``, and ``vibrosim_demo3_comsol.m``. The COMSOL file is a
matlab script to set up a COMSOL model for use in Vibrothermographic testing.
Many ``processtrak`` steps in VibroSim Simulator will interact with such a
model, however VibroSim Simulator itself is platform agnostic. The description
of the script itself will be left for the sister package ``VibroSim COMSOL``.

``vibrosim_demo3.prx``
---------------------

The ``.prx`` file contains the processing steps to be performed with the model. It is managed by the software tool called ``processtrak``, a part of ``Limatix``.

Use the tool in the following way. 

``processtrak <args> vibrosim_demo3.prx``

The following are possible arguments.

1. ``-s <steps>`` : Run only listed steps (multiple OK) 
2. ``-a`` : Run all steps
3. ``-i`` : Use ipython interactive mode to execute script
4. ``--needed`` : Filter down steps and input files according to what "needs"
    to be run -- i.e. missing or out-of-order steps, etc. DOES NOT PERFORM
    PROVENANCE VERIFICATION

Running this command: 

``processtrak --status vibrosim_demo3.prx``

produces the following output::

    Input file: vibrosim_demo3.xlg
    ---------------------------
               copyinput NOT_EXECUTED NEEDED
             dummyoutput NOT_EXECUTED NEEDED
              buildmodel NOT_EXECUTED NEEDED
                runmodal NOT_EXECUTED NEEDED
      synthetic_spectrum NOT_EXECUTED NEEDED
         entersweepfreqs NOT_EXECUTED NEEDED
     setsweepfrequencies NOT_EXECUTED NEEDED
          sweep_analysis NOT_EXECUTED NEEDED
          enterburstfreq NOT_EXECUTED NEEDED
       setburstfrequency NOT_EXECUTED NEEDED
          burst_analysis NOT_EXECUTED NEEDED
    eval_closure_state_from_tip_positions NOT_EXECUTED NEEDED
    calc_heating_singlefrequency NOT_EXECUTED NEEDED
       heatflow_analysis NOT_EXECUTED NEEDED

None of the steps have been run yet, so they all have the same ``NOT_EXECUTED`` and ``NEEDED`` tags. 

Now we can run the first step with the following command:

``processtrak vibrosim_demo3.prx -s copyinput``::

    Processing step copyinput on vibrosim_demo3.xlg->vibrosim_demo3.xlp

``processtrak --status vibrosim_demo3.prx``::

	Input file: vibrosim_demo3.xlg
	---------------------------
	           copyinput     2020-06-18T13:48:49.707642-05:00 
	         dummyoutput NOT_EXECUTED NEEDED
	          buildmodel NOT_EXECUTED NEEDED
	            runmodal NOT_EXECUTED NEEDED
	  synthetic_spectrum NOT_EXECUTED NEEDED
	     entersweepfreqs NOT_EXECUTED NEEDED
	 setsweepfrequencies NOT_EXECUTED NEEDED
	      sweep_analysis NOT_EXECUTED NEEDED
	      enterburstfreq NOT_EXECUTED NEEDED
	   setburstfrequency NOT_EXECUTED NEEDED
	      burst_analysis NOT_EXECUTED NEEDED
	eval_closure_state_from_tip_positions NOT_EXECUTED NEEDED
	calc_heating_singlefrequency NOT_EXECUTED NEEDED
	   heatflow_analysis NOT_EXECUTED NEEDED

Now that the ``copyinput`` step has been executed, it is no longer ``needed``
and the timestamp of the step has been recorded. This data is kept in the
``.xlp`` file.

Steps can be run out of order, as long as the ``.xlp`` has everything that is needed for the step.

``processtrak vibrosim_demo3.prx -s entersweepfreqs``

``processtrak --status vibrosim_demo3.prx``::

	Input file: vibrosim_demo3.xlg
	---------------------------
	           copyinput     2020-06-18T13:48:49.707642-05:00 
	         dummyoutput NOT_EXECUTED NEEDED
	          buildmodel NOT_EXECUTED NEEDED
	            runmodal NOT_EXECUTED NEEDED
	  synthetic_spectrum NOT_EXECUTED NEEDED
         entersweepfreqs     2020-06-18T13:53:03.368990-05:00  FAILURE NEEDED
	 setsweepfrequencies NOT_EXECUTED NEEDED
	      sweep_analysis NOT_EXECUTED NEEDED
	      enterburstfreq NOT_EXECUTED NEEDED
	   setburstfrequency NOT_EXECUTED NEEDED
	      burst_analysis NOT_EXECUTED NEEDED
	eval_closure_state_from_tip_positions NOT_EXECUTED NEEDED
	calc_heating_singlefrequency NOT_EXECUTED NEEDED
	   heatflow_analysis NOT_EXECUTED NEEDED

All steps will be run from scratch if the following command is run:

``processtrak vibrosim_demo3.prx -a``

All steps with the ``NEEDED`` flag will be run with the following command:

``processtrak vibrosim_demo3.prx -a --needed``


``vibrosim_demo3.xlg``
---------------------

A ``.xlg`` contains the unprocessed experiment log. ``processtrak`` is
primarily a tool for processing data collected in an experiment, after all. It
contains all the parameters necessary to start the simulation. The first step
in a simulation is to copy the ``.xlg`` into a processed experiment log
``.xlp`` file. This new ``.xlp`` file contains all the parameters in the
``.xlg`` and can be added to with ``processtrak`` steps. 

``.xlg`` files are ``xml`` based, meaning they are hierarchical in nature. It
is an single ``experiment`` tag with multiple ``measurement`` tags. Typically
there is a single ``measurement``, and simulations with different parameters
should be split into different ``.prx`` and ``.xlg`` files.

``submeas`` tags can be used to run multiple simulations at the same time with
mostly the same parameters save a few that are included in the ``submeas`` tag.
For example::

    <dc:experiment xmlns:dc="http://limatix.org/datacollect" xmlns:xlink="http://www.w3.org/1999/xlink" xmlns:dcv="http://limatix.org/dcvalue" xmlns:prx="http://limatix.org/processtrak/processinginstructions">
        <dc:measurement>
            <dc:measident>meas1</dc:measident> <!-- measident is used as a filename prefix for the various output files generated -->
            ...
            <dc:submeas>
              <dc:measident>meas1_direct_singlefreq</dc:measident>
              <dc:heatcalctype>singlefrequency</dc:heatcalctype>
            </dc:submeas>
            <dc:submeas>
              <dc:measident>meas1_via_weldercalc</dc:measident>
              <dc:heatcalctype>welder</dc:heatcalctype>
            </dc:submeas>
        </dc:measurement>
    </dc:experiment>  

In this example ``dc:measident`` is defined in both the ``measurement`` and the
``submeas`` levels. ``processtrak`` steps that are instructed to look for
submeasurements will pay attention to these changes.

This ``processtrak`` step from ``vibrosim_demo3_crosscheck.prx`` will only act on
submeasurements where the ``dc:heatcalctype`` is equal to
``'singlefrequency'``.::

  <prx:step name="calc_heating_singlefrequency">
    <prx:elementmatch>dc:measurement/dc:submeas[dc:heatcalctype='singlefrequency']</prx:elementmatch>
    <prx:script name="vibrosim_calc_heating_singlefrequency.py"/>
  </prx:step>

Alternatively, this ``processtrak`` step will act on all submeasurements.::
    
    <prx:step name="heatflow_analysis">
      <prx:elementmatch>dc:measurement/dc:submeas</prx:elementmatch>
      <prx:script name="vibrocomsol_heatflow_analysis_comsol.m"/>
    </prx:step>
