%This code will create a VibroSim model for an angled specimen using a completely
%imported geometry (as a .step file)

%This should be the first line in a VibroSim script. 
[M,model]=InitializeVibroSimScript();




%%%%%List necessary parameters
%%call the Vibrosim default parameters
VibroSim_default_params(M)
%% Setting up specimen geometries C14-UTCA-016W
%variables for specimen and isolators
laserfromedge_length = 0.0785 %m
isolator_length = 0.038; % m
isolator_width = 0.003; %
%it is clamped on
%LeftIsolator = 0.035; %Will be the location that it is clamped at.(~antinode for 5th mode)
%RightIsolator = spclength - 0.035;
isolator_dashpot_per_area = 670000 ;
% Change the excitation location to quarter from the side, excite torsional modes in the simulation.
couplant_coord=[0.025, 0.0048, 0.078, NaN    ]; 
% Define isolator locations, this is the center of the isolator
isolator_coords = [ 0.0319, 0, 0.0265 , 0.0  ,  % Bottom 
                    0.0319, 0.0048, 0.0265 , 0.0  ]; % Top


%%%%%Add parameters to the Param database
% Laser (displacement or velocity detection) coordinates
AddParamToParamdb(M,'lasery',0.0048,'m');
AddParamToParamdb(M,'laserx',0.025,'m');
AddParamToParamdb(M,'laserz',0.100,'m');
AddParamToParamdb(M,'laserdy',1);
AddParamToParamdb(M,'laserdx',0);
AddParamToParamdb(M,'laserdz',0);
% Specimen Parameters
AddParamToParamdb(M,'spcmaterialdampingtype','other');

use_impulse_force_excitation=true

%center of eigenvalues search
AddParamToParamdb(M,'simulationeigsaround',10000,'Hz')
% Freq Domain Study harmonicsweep Parameters
AddParamToParamdb(M,'simulationfreqstart',7000,'Hz');
AddParamToParamdb(M,'simulationfreqstep',5,'Hz');
AddParamToParamdb(M,'simulationfreqend',10000,'Hz');
AddParamToParamdb(M,'xducerforce',0.0,'N');
AddParamToParamdb(M,'staticload_mount',0.0,'N');
AddParamToParamdb(M,'xducerdisplacement',0.01,'m');
ObtainDCParameter(M,'xducerdisplacement','m');

% Isolator Parameters 
AddParamToParamdb(M,'isolatordashpotcoeff',isolator_dashpot_per_area,'Pa*s/m'); %Ns/m * 1/m^2
AddParamToParamdb(M,'isolatorthickness', 0.001, 'm');     % Z Direction
AddParamToParamdb(M,'isolatorlength'   , isolator_length, 'm');       % X Direction
AddParamToParamdb(M,'isolatorwidth'    , isolator_width,'m');     % Y Direction
% Couplant parameters
AddParamToParamdb(M,'couplantlength',.003,'m');
AddParamToParamdb(M,'couplantwidth',.003,'m');
AddParamToParamdb(M,'couplantdashpotcoeff',0,'Pa*s/m'); %Ns/m * 1/m^2
% Time domain study timings
AddParamToParamdb(M,'timedomain_start_time',-4e-6,'s');
AddParamToParamdb(M,'timedomain_step_time',1e-6,'s');
AddParamToParamdb(M,'timedomain_end_time',20e-3,'s'); % Double the time so I can fit sinusoids to the data
% Mesh study parameters
meshsizemin = 1e-3 ;
meshsize = 4e-3 ; 
AddParamToParamdb(M,'spcmeshsize'     ,meshsize,'m');
AddParamToParamdb(M,'isolatormeshsize',meshsize,'m');
AddParamToParamdb(M,'couplantmeshsize',meshsize,'m');


%%%%%Define and build the geometry
% Define a procedure for building the geometry. Steps can be sequenced by using
% the pipe (vertical bar | ) character. 
bldgeom = @(M,geom) angle_geometry(M,geom,'specimen',meshsizemin,meshsize) | ...
	  @(specimen) AttachThinCouplantIsolators(M,geom,specimen,couplant_coord,isolator_coords); 
	  

% Define a procedure for creating the various physics definitions. Steps can be 
% sequenced by using the pipe (vertical bar | ) character. 
bldphysics = @(M,geom,specimen,flaw) VibroPhysics(M,geom,specimen,flaw,'modal|harmonicsweep');

% Define a procedure to create needed COMSOL result nodes. Generally can pass
% either the wrapped model (M) or unwrapped (model). Pipelining will only work
% with the wrapped model 'M'.
genresults = @(model) VibroResults(model),

% Define a path and filename for saving the generated COMSOL model.
% BuildVibroModel() will create the model's directory if necessary
filename = 'dashpot_impulse_response.mph'; 
full_path = mfilename('fullpath');
path_bits = strsplit(full_path,'/');
path_part = path_bits(1:end-1);

disp('file path')
savefilename  = fullfile('/',path_part{:},filename)

% Given the model wrapper, the procedures for building the geometry, the crack, the physics,
% and for storing the results, and a file name for saving, have COMSOL build the model.x
BuildVibroModel(M,...
		bldgeom, ...
		[], ...
		bldphysics, ...
		genresults, ...
		savefilename);

RunAllStudies(model);
mphsave(M.node,savefilename);		
