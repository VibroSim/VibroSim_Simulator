% function ret = vibrosim_demo3_comsol(dc_dest_href,dc_measident_str,dc_dummy_heatingdata_href,dc_amplitude_float,dc_staticload_mount_float,dc_spcmaterial_str,dc_YoungsModulus_float, dc_YieldStrength_float, dc_PoissonsRatio_float, dc_Density_float,dc_spcThermalConductivity_float, dc_spcSpecificHeatCapacity_float,dc_spcrayleighdamping_alpha_float,dc_spcrayleighdamping_beta_float, dc_exc_t0_float, dc_exc_t4_float, dc_simulationcameranetd_float,dc_cracksemimajoraxislen_float,dc_cracksemiminoraxislen_float,dc_crack_type_side1_str,dc_crack_type_side2_str)
%> @brief Here is a third example of how to use BuildVibroModel
%> The commented function declaration is used by processtrak to
%> figure out the parameters to pass.

% NOTE: The first line (commented out function definition)
%       is important because it tells processtrak which parameters
%       to pass
%
% NOTE: The name of this script should end with "_comsol.m" as that
%       makes processtrak run it with COMSOL MATLAB rather than
%       bare MATLAB>


% This is identical to the second example, except static mount boundary conditions are added
% using 'Seq' and it automatically runs as well as builds
%
% NOTE: Make your window as wide as possible when viewing this file!

%                 x          y         z      angle
couplant_coord=[ .14/2,   .0254/2,     .012,      NaN    ];
isolator_coords=[.13,     .0254/2,     0,   0.0,  % top-left
		 .01,     .0254/2,     0,   0.0,  % top-right
                 .12,     .0254/2,     .012,      0.0,  % bottom-left
		 .02,     .0254/2,     .012,      0.0]; % bottom-right

% ( NaN for angle causes it to create a circular isolator)

% You may wish to uncomment these next two lines if this is part of a function
% -- that way if the function fails you can access the wrapped and unwrapped model 
% variables just because they are globals. 
%global M
%global model


% InitializeVibroSimScript() connects to COMSOL and initializes and 
% returns the wrapped model variable (M) and the unwrapped node (model). 
[M,model]=InitializeVibroSimScript();


 
% Call a function that sets various parameters to be used by the model. 
% Defaults specified in this function can be overridden by values set
% below
VibroSim_default_params(M);



% Set COMSOL parameters from experiment log, passed as parameters
AddParamToParamdb(M,'spcYoungsModulus',dc_YoungsModulus_float,'Pa');
AddParamToParamdb(M,'spcPoissonsRatio',dc_PoissonsRatio_float,'');
AddParamToParamdb(M,'spcDensity',dc_Density_float,'kg/m^3');
AddParamToParamdb(M,'spcThermalConductivity',dc_spcThermalConductivity_float,'W/m/K');
AddParamToParamdb(M,'spcSpecificHeatCapacity',dc_spcSpecificHeatCapacity_float,'J/kg/K');

% simulationtimestart, simulationtimestep, and simulationtimeend specify the time range of the heat flow simulation
AddParamToParamdb(M,'simulationtimestart',dc_exc_t0_float,'s');
AddParamToParamdb(M,'simulationtimestep',0.02,'s');
AddParamToParamdb(M,'simulationtimeend',dc_exc_t4_float+0.8,'s'); % .8 seconds after assumed end of vibration

AddParamToParamdb(M,'spcmaterialdampingtype','RayleighDamping');
AddParamToParamdb(M,'spcrayleighdamping_alpha',dc_spcrayleighdamping_alpha_float,'s^-1');
AddParamToParamdb(M,'spcrayleighdamping_beta',dc_spcrayleighdamping_beta_float,'s');




% Crack position
crackx = '.07[m]';
cracky = '.0127[m]';
crackz = '0[m]';



% Laser (displacement or velocity detection) coordinates
AddParamToParamdb(M,'laserx',.07,'m');
AddParamToParamdb(M,'lasery',.0254/4.0,'m');
AddParamToParamdb(M,'laserz',0.0,'m');

AddParamToParamdb(M,'laserdx',0);
AddParamToParamdb(M,'laserdy',0);
AddParamToParamdb(M,'laserdz',1);


% Because we are modeling linear single-frequency excitation, we need to 
% select a transducer calibration file. 

% Load in xducercalib file to xducercalib function, set up xducerdisplacement as variable (WARNING: This step can be slow!)
CreateTransducerDisplacementVariable(M,fullfile(fileparts(which('BuildVibroModel')),'..','transducer_data','constant_10micronpervolt_displacementampl.dat'),dc_amplitude_float);


% Camera noise parameter
CreateCameraNoise(M,'cameranoise',dc_simulationcameranetd_float);




% Crack type --> 'penny' or 'through'
cracktype = 'penny';
assert(contains(dc_crack_type_side1_str,cracktype) | strcmpi(dc_crack_type_side1_str,'none')) % crack side types must be consistent with geometry being built 
assert(contains(dc_crack_type_side2_str,cracktype) | strcmpi(dc_crack_type_side2_str,'none')) % crack side types must be consistent with geometry being built 


% staticload_mount and xducerforce are not really needed because we're not doing the static analysis. 
% Nevertheless, the static analysis is configured so we still have to set them
AddParamToParamdb(M,'staticload_mount',0.0,'N');
AddParamToParamdb(M,'xducerforce',0.0,'N');
ObtainDCParameter(M,'staticload_mount','N');
ObtainDCParameter(M,'xducerforce','N');

% Define a procedure for building the geometry. Steps can be sequenced by using
% the pipe (vertical bar | ) character. 
  bldgeom = @(M,geom) CreateRectangularBarSpecimen(M,geom,'specimen', '.14[m]', '.0254[m]', '.012[m]',dc_spcmaterial_str) | ... % length, width, thickness, material
	  @(specimen) AttachThinCouplantIsolators(M,geom,specimen, ...
						    couplant_coord, ...
						    isolator_coords) | ...
	  ... % Top-left and top-right isolators (.isolators{1} and .isolators{2}: net force of -staticload_mount in the z direction
	  @(specimen) AddBoundaryCondition(M,specimen,specimen.isolators{1},...
					   [specimen.isolators{1}.tag '_tltrforce'], ...
					   'solidmech_static',...  % physics
					   'staticloading', ...    % BC class
					   @(M,physics,bcobj) ...
					    BuildFaceTotalForceBC(M,geom,physics,specimen.isolators{1},bcobj,...
								  @(M,geom,tlisolator) ...
								   union(specimen.isolators{1}.getfreefaceselection(M,geom,specimen.isolators{1}), ...
									 specimen.isolators{2}.getfreefaceselection(M,geom,specimen.isolators{2})), ...
								  { 0, 0 , -dc_staticload_mount_float })) | ...
          ... % Bottom-left isolator (.isolators{3}) is fixed
          @(specimen) AddBoundaryCondition(M,specimen,specimen.isolators{3}, ...
					   [specimen.isolators{3}.tag '_blfixed'],...
					   'solidmech_static', ...  % physics
					   'staticloading', ...     % BC class
					   @(M,physics,bcobj) ...
					    BuildFaceFixedBC(M,geom,physics, ...
							     specimen.isolators{3}, ...
							     bcobj, ...
							     specimen.isolators{3}.getfreefaceselection)) |  ...
	  ... % Bottom-right isolator (.isolators{4}) is rolling -- displacement in the Z direction is zero
	  @(specimen) AddBoundaryCondition(M,specimen,specimen.isolators{4},[specimen.isolators{4}.tag '_brrolling'], ...
					   'solidmech_static', ...
					   'staticloading', ...
					   @(M,physics,bcobj) ...
					    BuildFaceDirectionalDisplacementBC(M,geom,physics,...
									       specimen.isolators{4},...
									       bcobj, ...
									       specimen.isolators{4}.getfreefaceselection, ...
									       GetOutwardNormal(M,geom,specimen.isolators{4}.pos,specimen.isolators{4}.centerpos), ...
									       0.0)) | ...
	  ... % Force condition on couplant
	  @(specimen) AddBoundaryCondition(M,specimen,specimen.couplant,[specimen.couplant.tag '_xducerforce'], ...
					   'solidmech_static', ...
					   'staticloading', ...
					   @(M,physics,bcobj) ...
					    BuildFaceTotalForceBC(M,geom,physics,specimen.couplant,bcobj, ...
								  specimen.couplant.getfreefaceselection, ...
								  ... % xducerforce is multiplied by outward normal 
		  			                          ... % on the couplant, which is equivalent to the
					                          ... % inward normal on the specimen
								  MultiplyScalarStrByNumericVec('xducerforce', ...
												GetOutwardNormal(M,geom,specimen.couplant.pos,specimen.couplant.centerpos)))) | ...
	  @(specimen) AddView(M,specimen,11.54, ... % zoom angle
			      { '-.32[m]','-.5[m]','-.39[m]' }, ...% camera position
			      { crackx, ...  % camera target -- same location as crack
				cracky, ...
				crackz }, ...
			      { '.3087','.4116','-.857' }, ... % up direction
			      { crackx, ... % camera rotation point
				cracky, ...
				crackz });



% Define a procedure for building the crack. Steps can be sequenced by using
% the pipe (vertical bar | ) character. 
bldcrack = @(M,geom,specimen) CreateCrack(M,geom,'crack',specimen, ...
					  { crackx, ... % crack location
					    cracky, ...
					    crackz }, ...
					  dc_cracksemimajoraxislen_float, ...
					  dc_cracksemiminoraxislen_float, ...
					  [0,1,0], ... % Axis major direction (surface growth)
					  [0,0,-1], ... % Axis minor direction (depth growth for a surface crack)
					  [ .001, .002, .003 ], ... % Subradii
                                          {'solidmech_harmonicsweep','solidmech_harmonicburst'}, ...
					  dc_dummy_heatingdata_href{1},... % Text file to hold crack heating energies. 
                      cracktype); 

% Define a procedure for creating the various physics definitions. Steps can be 
% sequenced by using the pipe (vertical bar | ) character. 
bldphysics = @(M,geom,specimen,flaw) VibroPhysics(M,geom,specimen,flaw,'broadbandprocess');

% Define a procedure to create needed COMSOL result nodes. Generally can pass
% either the wrapped model (M) or unwrapped (model). Pipelining will only work
% with the wrapped model 'M'.
genresults = @(M) VibroResults(M);

% Define a path and filename for saving the generated COMSOL model.
% BuildVibroModel() will create the model's directory if necessary
%savefilename = fullfile(tempdir,sprintf('vibrosim_%s',char(java.lang.System.getProperty('user.name'))),'vibrosim_demo3.mph');

savefilename = fullfile(dc_dest_href{1},[ dc_measident_str '.mph']);


% Given the model wrapper, the procedures for building the geometry, the crack, the physics,
% and for storing the results, and a file name for saving, have COMSOL build the model.x
BuildVibroModel(M,...
		bldgeom, ...
		bldcrack, ...
		bldphysics, ...
		genresults, ...
		savefilename);

% Uncomment these next two lines to run the model and save the output after building the model
%RunAllStudies(model);
%mphsave(M.node,savefilename);

ret = { { 'dc:model_comsol', { savefilename } } };
