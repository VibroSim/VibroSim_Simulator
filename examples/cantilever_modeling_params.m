%> You can define parameters using AddParamToParamdb.
%> As an alternative you can define COMSOL parameters with CreateParameter()
%> As another alternative you can store parameters in the M structure
%> (fields can be added to M, with addprop()).

function default_params(M)

% crack parameters
AddParamToParamdb(M,'cracksemimajoraxislen',3e-3,'m');
AddParamToParamdb(M,'cracksemiminoraxislen',1.5e-3,'m');


% multisweep parameters: 11/9... These are overridden by the set_freqbands step
AddParamToParamdb(M,'seg1_freqstart',0.000,'Hz');
AddParamToParamdb(M,'seg1_freqstep',0.476837158203125,'Hz');
AddParamToParamdb(M,'seg1_freqend',74.3865966796875,'Hz');

AddParamToParamdb(M,'seg2_freqstart',73.909759521484375,'Hz');
AddParamToParamdb(M,'seg2_freqstep',7.62939453125,'Hz');
AddParamToParamdb(M,'seg2_freqend',2660.27450561523438,'Hz');

AddParamToParamdb(M,'seg3_freqstart',2662.65869140625,'Hz');
AddParamToParamdb(M,'seg3_freqstep',61.03515625,'Hz');
AddParamToParamdb(M,'seg3_freqend',23231.50634765625,'Hz');

AddParamToParamdb(M,'seg4_freqstart',23236.7515563964844,'Hz');
AddParamToParamdb(M,'seg4_freqstep',488.28125,'Hz');
AddParamToParamdb(M,'seg4_freqend',999999.523162841797,'Hz');



%AddParamToParamdb(M,'staticload_mount',10000,'N');
%AddParamToParamdb(M,'xducerforce',100,'N');
AddParamToParamdb(M,'simulationfreqstart',10,'Hz');
AddParamToParamdb(M,'simulationfreqstep',10,'Hz');
AddParamToParamdb(M,'simulationfreqend',90000,'Hz');
AddParamToParamdb(M,'simulationburstfreq',15270,'Hz');
AddParamToParamdb(M,'simulationneigs',70);  % Number of frequencies for modal analysis to seek out
AddParamToParamdb(M,'simulationeigsaround',100,'Hz'); % Center frequency for modal analysis

AddParamToParamdb(M,'simulationtimestart',0.2,'s');
AddParamToParamdb(M,'simulationtimestep',0.02,'s');
AddParamToParamdb(M,'simulationtimeend',1.8,'s');

AddParamToParamdb(M,'simulationcameranetd',.022,'K');
AddParamToParamdb(M,'simulationsurfaceemissivity',1.0);
ObtainDCParameter(M,'simulationsurfaceemissivity');

% Laser (displacement or velocity detection) coordinates
AddParamToParamdb(M,'laserx',.07,'m');
AddParamToParamdb(M,'lasery',.0254/4.0,'m');
AddParamToParamdb(M,'laserz',0.0,'m');

% Laser (displacement or velocity detection) direction vector
AddParamToParamdb(M,'laserdx',0);
AddParamToParamdb(M,'laserdy',0);
AddParamToParamdb(M,'laserdz',1);

% Crack position
AddParamToParamdb(M,'simulationcrackx',.07,'m');
AddParamToParamdb(M,'simulationcracky',.0254/2.0,'m');
AddParamToParamdb(M,'simulationcrackz',0,'m');

% Crack heating model parameters
%AddParamToParamdb(M,'expV0',20000e6,'W*Pa/m^2/Hz'); % Power law strain dependence: leading coefficient
%AddParamToParamdb(M,'V1',1.7); % Power law strain dependence: exponent
%AddParamToParamdb(M,'m1',0.00152,'m')
%AddParamToParamdb(M,'l0',1.0/0.000762,'1/m')
AddParamToParamdb(M,'wh',20e6,'Pa');

% Crack heating model requires these parameters instantiated
%ObtainDCParameter(M,'expV0','W*Pa/m^2/Hz');
%ObtainDCParameter(M,'V1');
%ObtainDCParameter(M,'m1','m');
%ObtainDCParameter(M,'l0','1/m');
ObtainDCParameter(M,'wh','Pa');


AddParamToParamdb(M,'amplitude',1.0,'V');
AddParamToParamdb(M,'spclength',0.2538,'m'); % length, measured by Jared Taylor
AddParamToParamdb(M,'spcwidth',50.9e-3,'m'); % width, measured by Jared Taylor
AddParamToParamdb(M,'spcthickness',6.45e-3,'m'); % thickness, measured by Jared Taylor


% These parameters are book values for Aluminum (6061?) material 
%AddParamToParamdb(M,'spcmaterial','Aluminum');
%AddParamToParamdb(M,'spcYoungsModulus',69e9,'Pa');
%AddParamToParamdb(M,'spcDensity',2697,'kg/m^3'); % Density, measured by Jared Taylor
%AddParamToParamdb(M,'spcPoissonsRatio',0.32,''); % book value
AddParamToParamdb(M,'spcThermalConductivity',167,'W/m/K');
AddParamToParamdb(M,'spcSpecificHeatCapacity',896,'J/kg/K');

%% These parameters are book values for Inconel 718  material 
%AddParamToParamdb(M,'spcmaterial','Inconel-718');
%AddParamToParamdb(M,'spcYoungsModulus',200e9,'Pa');
%AddParamToParamdb(M,'spcDensity',8190,'kg/m^3');
%AddParamToParamdb(M,'spcPoissonsRatio',0.33,''); % book value
%AddParamToParamdb(M,'spcEta',1e-5,'');
%AddParamToParamdb(M,'spcThermalConductivity',11.4,'W/m/K');
%AddParamToParamdb(M,'spcSpecificHeatCapacity',435.0,'J/kg/K');

AddParamToParamdb(M,'spcmeshtype','TETRAHEDRAL');
%AddParamToParamdb(M,'spcmeshtype','HEXAHEDRAL');
AddParamToParamdb(M,'spcmeshsize',.004,'m');
%AddParamToParamdb(M,'spcfacemethod','FreeQuad');
%AddParamToParamdb(M,'spcsweepelements',15)

AddParamToParamdb(M,'tlmountoffsetx',.13,'m');
AddParamToParamdb(M,'blmountoffsetx',.12,'m');
AddParamToParamdb(M,'brmountoffsetx',.02,'m');
AddParamToParamdb(M,'trmountoffsetx',.01,'m');
AddParamToParamdb(M,'xduceroffsetx',.07,'m');



% NOTE: isolatornlayers removed -- assume something else 
% is doing the multiplication to determine total isolatorthickness
AddParamToParamdb(M,'isolatorthickness',.00125,'m');
AddParamToParamdb(M,'isolatorlength',.010,'m');
AddParamToParamdb(M,'isolatorwidth',.0254,'m');
AddParamToParamdb(M,'isolatormeshtype','HEXAHEDRAL');
AddParamToParamdb(M,'isolatorfacemethod','FreeQuad');
AddParamToParamdb(M,'isolatormeshsize',.004,'m');
AddParamToParamdb(M,'isolatorsweepelements',6,'');
AddParamToParamdb(M,'isolatorThermalConductivity',.05,'W/m/K');  % this number may not be very meaningful
AddParamToParamdb(M,'isolatorSpecificHeatCapacity',2500,'J/kg/K');  % this number may not be very meaningful

% WARNING: These values do not correspond to any particular material
AddParamToParamdb(M,'isolatormaterial','isocardstock');
AddParamToParamdb(M,'isolatorYoungsModulus',3.5e6,'Pa');
AddParamToParamdb(M,'isolatorDensity',870,'kg/m^3');
AddParamToParamdb(M,'isolatorPoissonsRatio',0.3,'');
% 6/11/18 change isolator damping behavior from Eta to d_A.... where d_A=(spring_constant_per_unit_area/(2*pi*f))*eta  
% This works in time domain and is more physically consistent with observations (behaves like dashpot)
% !!!*** Should do the same for couplant and perhaps specimen
  % k/area = isolatorYoungsModulus/isolatorthickness = 3.5e6/.00125 = 2.8e9 (N/m)/m^2.. at a given frequency
  % equivalent dashpot coeff. for equivalent impedance magnitude is k/omega... at 1 kHz = 2.8e9/(2*pi*1e3)=4.5e5 (N*s/m)/m^2
%AddParamToParamdb(M,'isolatorEta',.5,'');
%AddParamToParamdb(M,'isolatordashpotcoeff',11e3,'Pa*s/m');
% 1/30/19 change isolatordashpotcoeff from 11e3 Pa*s/m to 100e3
AddParamToParamdb(M,'isolatordashpotcoeff',100e3,'Pa*s/m');

AddParamToParamdb(M,'couplantthickness',.00025,'m');
AddParamToParamdb(M,'couplantlength',.010,'m');
AddParamToParamdb(M,'couplantwidth',.010,'m');
AddParamToParamdb(M,'couplantmeshtype','HEXAHEDRAL');
AddParamToParamdb(M,'couplantfacemethod','FreeQuad');
AddParamToParamdb(M,'couplantmeshsize',.004,'m');
AddParamToParamdb(M,'couplantsweepelements',4,'');

% WARNING: These values do not correspond to any particular material
AddParamToParamdb(M,'couplantmaterial','coupcardstock');
%AddParamToParamdb(M,'couplantYoungsModulus',3.5e6,'Pa');
AddParamToParamdb(M,'couplantYoungsModulus',0,'Pa');
AddParamToParamdb(M,'couplantDensity',870,'kg/m^3');
AddParamToParamdb(M,'couplantPoissonsRatio',0.3,'');
%AddParamToParamdb(M,'couplantdashpotcoeff',100e3,'Pa*s/m');
AddParamToParamdb(M,'couplantdashpotcoeff',0.0,'Pa*s/m');
AddParamToParamdb(M,'couplantThermalConductivity',.05,'W/m/K');  % this number may not be very meaningful
AddParamToParamdb(M,'couplantSpecificHeatCapacity',2500,'J/kg/K');  % this number may not be very meaningful


AddParamToParamdb(M,'meshsizemin',.001,'m');
AddParamToParamdb(M,'meshsize',.004,'m');

AddParamToParamdb(M,'excitation_t0',0.2,'s');
AddParamToParamdb(M,'excitation_t1',0.21,'s');
AddParamToParamdb(M,'excitation_t2',1.199,'s');
AddParamToParamdb(M,'excitation_t3',1.2,'s');

AddParamToParamdb(M,'impulseexcitation_t0',0.0,'s');
AddParamToParamdb(M,'impulseexcitation_width',2.0e-6,'s');

AddParamToParamdb(M,'timedomain_start_time',-4e-6,'s');
AddParamToParamdb(M,'timedomain_step_time',1e-6,'s');
AddParamToParamdb(M,'timedomain_end_time',50e-3,'s');


CreateExcitationWindow(M,'excitationwindow', ...
		       ObtainDCParameter(M,'excitation_t0','s'), ...
		       ObtainDCParameter(M,'excitation_t1','s'), ...
		       ObtainDCParameter(M,'excitation_t2','s'), ...
		       ObtainDCParameter(M,'excitation_t3','s'));

CreateImpulseExcitation(M,'impulse_excitation', ...
			ObtainDCParameter(M,'impulseexcitation_t0','s'), ...
			ObtainDCParameter(M,'impulseexcitation_width','s'));

CreateCameraNoise(M,'cameranoise',ObtainDCParameter(M,'simulationcameranetd','K'));

