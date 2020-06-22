%> You can define parameters using AddParamToParamdb.
%> As an alternative you can define COMSOL parameters with CreateParameter()
%> As another alternative you can store parameters in the M structure
%> (fields can be added to M, with addprop()).

function cantilever_modeling_params(M)



% Dummy segment boundaries -- generally set by vibrocomsol_set_freqbands_comsol
AddParamToParamdb(M,'seg1_freqstart',0.000,'Hz');
AddParamToParamdb(M,'seg1_freqstep',1.0,'Hz');
AddParamToParamdb(M,'seg1_freqend',0.0,'Hz');

AddParamToParamdb(M,'seg2_freqstart',0.0,'Hz');
AddParamToParamdb(M,'seg2_freqstep',1.0,'Hz');
AddParamToParamdb(M,'seg2_freqend',0.0,'Hz');

AddParamToParamdb(M,'seg3_freqstart',0.0,'Hz');
AddParamToParamdb(M,'seg3_freqstep',1.0,'Hz');
AddParamToParamdb(M,'seg3_freqend',0.0,'Hz');

AddParamToParamdb(M,'seg4_freqstart',0.0,'Hz');
AddParamToParamdb(M,'seg4_freqstep',1.0,'Hz');
AddParamToParamdb(M,'seg4_freqend',0.0,'Hz');


AddParamToParamdb(M,'simulationfreqstart',10,'Hz');
AddParamToParamdb(M,'simulationfreqstep',10,'Hz');
AddParamToParamdb(M,'simulationfreqend',90000,'Hz');
AddParamToParamdb(M,'simulationburstfreq',15270,'Hz');
AddParamToParamdb(M,'simulationneigs',70);  % Number of frequencies for modal analysis to seek out
AddParamToParamdb(M,'simulationeigsaround',100,'Hz'); % Center frequency for modal analysis


AddParamToParamdb(M,'simulationsurfaceemissivity',1.0);
ObtainDCParameter(M,'simulationsurfaceemissivity');



AddParamToParamdb(M,'spcmeshtype','TETRAHEDRAL');
%AddParamToParamdb(M,'spcmeshtype','HEXAHEDRAL');
AddParamToParamdb(M,'spcmeshsize',.004,'m');
%AddParamToParamdb(M,'spcfacemethod','FreeQuad');
%AddParamToParamdb(M,'spcsweepelements',15)





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


