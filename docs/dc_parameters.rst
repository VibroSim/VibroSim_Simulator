DC (DataCollect) Parameters
===========================

DC parameters are the parameters that are required for every vibrothermographic
simulation. They are defined in the ``.xlg`` and ``.xlp`` files. The examples
provided here are from the ``vibrosim_demo3`` example.

dest
-------
Where to save the generated output. ::

    <dc:dest xlink:type="simple" xlink:href="vibrosim_demo3_output/"/>

measident
---------
Used as a filename prefix for the various output files generated. It should be
redefined in submeasurements if they are being used.::

    <dc:measident>meas1</dc:measident> 

amplitude
---------
COULDN'T FIND IN USE, NOT IN THE \*_comsol.m file and not in the m_files directory.
text.::

    <dc:amplitude dcv:units="Volts">3</dc:amplitude>

friction_coefficient
--------------------
text. ::

    <dc:friction_coefficient>.0661</dc:friction_coefficient>

msqrtR
------
text. ::

    <dc:msqrtR dcv:units="m^-1.5">3.14e6</dc:msqrtR>

crack_model_shear_factor
------------------------
text. ::

    <dc:crack_model_shear_factor>1.0</dc:crack_model_shear_factor> 

staticload_mount
----------------
text. ::

    <dc:staticload_mount dcv:units="N">0</dc:staticload_mount> 

staticload
----------
text. ::

    <dc:staticload dcv:units="Pascals">0</dc:staticload> 

xducerforce
-----------
text. ::

    <dc:xducerforce dcv:units="Newtons">0</dc:xducerforce>

exclength
---------
text. ::

    <dc:exclength dcv:units="seconds">1.0</dc:exclength> <!-- Should eliminate this parameter... just use dc:exc_tx, below -->

tortuosity
----------
text. ::

    <dc:tortuosity dcv:units="degrees">30.5</dc:tortuosity><!-- 30.5 was the average from our training specimens -->

numdraws
--------
text. ::

    <dc:numdraws>100</dc:numdraws>

spcmaterial
-----------
text. ::

    <dc:spcmaterial>Ti 6-4</dc:spcmaterial>

YoungsModulus
-------------
text. ::

    <dc:YoungsModulus dcv:units="Pascals">117.9e9</dc:YoungsModulus>

YieldStrength
-------------
text. ::

    <dc:YieldStrength dcv:units="Pascals">944.58e6</dc:YieldStrength>

PoissonsRatio
-------------
text. ::

    <dc:PoissonsRatio>0.342</dc:PoissonsRatio>

Density
-------
text. ::

    <dc:Density dcv:units="kg/m^3">4430</dc:Density>

spcThermalConductivity
----------------------
text. ::

    <dc:spcThermalConductivity dcv:units="W/m/K">6.7</dc:spcThermalConductivity>

spcSpecificHeatCapacity
-----------------------
text. ::

    <dc:spcSpecificHeatCapacity dcv:units="J/kg/K">526.3</dc:spcSpecificHeatCapacity>

simulationcameranetd
--------------------
text. ::

    <dc:simulationcameranetd dcv:units="K">.022</dc:simulationcameranetd> <!-- magnitude of simulated camera noise: noise equivalent temperature difference (NETD) -->
    
spcrayleighdamping_alpha
------------------------
text. ::

    <dc:spcrayleighdamping_alpha dcv:units="s^-1">2</dc:spcrayleighdamping_alpha>

spcrayleighdamping_beta
-----------------------
text. ::

    <dc:spcrayleighdamping_beta dcv:units="s">5e-9</dc:spcrayleighdamping_beta>

spcviscousdamping
-----------------
text. ::

    <dc:spcviscousdamping dcv:units="N*s"> xxxx </dc:spcviscousdamping>
    
crack_model_normal
------------------
text. ::

    <dc:crack_model_normal>Tada_ModeI_CircularCrack_along_midline</dc:crack_model_normal>

crack_model_shear
-----------------
text. ::

    <dc:crack_model_shear>Fabrikant_ModeII_CircularCrack_along_midline</dc:crack_model_shear>

crack_shearstress_axis
----------------------
text. ::

    <dc:crack_shearstress_axis>major</dc:crack_shearstress_axis> 

crack_type_side1
----------------
text. ::

    <dc:crack_type_side1>quarterpenny</dc:crack_type_side1>

crack_type_side2
----------------
text. ::

    <dc:crack_type_side2>quarterpenny</dc:crack_type_side2>

thickness
---------
dc:thickness is used only for crack_type == halfthrough ::

    <dc:thickness dcv:units="meters">1e-3</dc:thickness>

reff_side1
----------
text. ::

    <dc:reff_side1 dcv:units="m" dcv:arraystorageorder="C"><dcv:arrayshape>9</dcv:arrayshape><dcv:arraydata>.5e-3 .7e-3 .9e-3 1.05e-3 1.2e-3 1.33e-3 1.45e-3 1.56e-3 1.66e-3</dcv:arraydata></dc:reff_side1>

seff_side1
----------
text. ::

    <dc:seff_side1 dcv:units="Pa" dcv:arraystorageorder="C"><dcv:arrayshape>9</dcv:arrayshape><dcv:arraydata>0.0 50e6 100e6 150e6 200e6 250e6 300e6 350e6 400e6</dcv:arraydata></dc:seff_side1>

reff_side2
----------
text. ::

    <dc:reff_side2 dcv:units="m" dcv:arraystorageorder="C"><dcv:arrayshape>9</dcv:arrayshape><dcv:arraydata> .5e-3 .7e-3 .9e-3 1.05e-3 1.2e-3 1.33e-3 1.45e-3 1.56e-3 1.66e-3</dcv:arraydata></dc:reff_side2>

seff_side2
----------
text. ::

    <dc:seff_side2 dcv:units="Pa" dcv:arraystorageorder="C"><dcv:arrayshape>9</dcv:arrayshape><dcv:arraydata> 0.0 50e6 100e6 150e6 200e6 250e6 300e6 350e6 400e6</dcv:arraydata></dc:seff_side2>

exc_t0
------
Start of envelope ramp-up. ::

    <dc:exc_t0 dcv:units="s">0.0</dc:exc_t0>  

exc_t1
------
End of enelope ramp-up. ::

    <dc:exc_t1 dcv:units="s">0.02</dc:exc_t1>

exc_t2
------
Start of envelope ramp-down ::

    <dc:exc_t2 dcv:units="s">0.98</dc:exc_t2>

exc_t3
------
End of envelope ramp-down. ::

    <dc:exc_t3 dcv:units="s">1.00</dc:exc_t3> 

exc_t4
------
End of vibration calculation. ::

    <dc:exc_t4 dcv:units="s">1.00</dc:exc_t4>
