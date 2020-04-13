% function ret = cantilever_example_comsol(dc_dest_href,dc_measident_str,dc_dummy_heatingdata_href,dc_xducerforce_float,dc_spcmaterial_str,dc_YoungsModulus_float, dc_YieldStrength_float, dc_PoissonsRatio_float, dc_Density_float,dc_spcThermalConductivity_float, dc_spcSpecificHeatCapacity_float,dc_spcrayleighdamping_alpha_float,dc_spcrayleighdamping_beta_float,dc_couplantx_float, dc_couplanty_float, dc_couplantz_float, dc_couplantangle_float, dc_exc_t0_float, dc_exc_t4_float, dc_simulationcameranetd_float)
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


%
% NOTE: Make your window as wide as possible when viewing this file!


% You may wish to uncomment these next two lines if this is part of a function
% -- that way if the function fails you can access the wrapped and unwrapped model 
% variables just because they are globals. 
%global M
%global model


% InitializeVibroSimScript() connects to COMSOL and initializes and 
% returns the wrapped model variable (M) and the unwrapped node (model). 
[M,model]=InitializeVibroSimScript();


 
% Call a function that sets various parameters to be used by the model. 
cantilever_modeling_params(M);


% Set COMSOL parameters from experiment log, passed as parameters
AddParamToParamdb(M,'xducerforce',dc_xducerforce_float,'N');
AddParamToParamdb(M,'spcmaterial',dc_spcmaterial_str);
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

% Camera noise parameter
CreateCameraNoise(M,'cameranoise',dc_simulationcameranetd_float);

%                 x          y         z      angle
%couplant_coord=[ .245,       .025,        0,      NaN    ];
couplant_coord=[ dc_couplantx_float,  dc_couplanty_float,  dc_couplantz_float,  dc_couplantangle_float   ];
isolator_coords=[];
% ( NaN for angle causes it to create a circular isolator)

  

% Extract needed parameters...

ObtainDCParameter(M,'xducerforce','N');



% Define a procedure for building the geometry. Steps can be sequenced by using
% the pipe (vertical bar | ) character. 
bldgeom = @(M,geom) CreateRectangularBarSpecimen(M,geom,'specimen') | ...
	  @(specimen) AttachThinCouplantIsolators(M,geom,specimen, ...
						    couplant_coord, ...
						    isolator_coords) | ...
	  ... % Add static boundary conditions to specimen
	  ... % Top-left and top-right isolators (.isolators{1} and .isolators{2}: net force of -staticload_mount in the z direction
	  @(specimen) AddBoundaryCondition(M,specimen,specimen,...
					   [specimen.tag '_fixedend'], ...
					   { 'solidmech_static','solidmech_harmonicper','solidmech_harmonic','solidmech_modal','solidmech_timedomain' },...  % physics
					   'fixed', ...    % BC class
					   @(M,physics,bcobj) ...
					    BuildFaceFixedBC(M,geom,physics,specimen,bcobj,...
							     @(M,geom,spec) GetBlockFace(M,geom,specimen,[-1,0,0]))) | ...
					    ... % BuildFaceSpringFoundationBC(M,geom,physics,specimen,bcobj,...
					    ... %				@(M,geom,spec) GetBlockFace(M,geom,specimen,[-1,0,0]), ...
					    ... %				... %FixedEnd_k_A, FixedEnd_DampPerArea)) |...
					    ... % 				FixedEnd_k_A_combined)) |...
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
												GetOutwardNormal(M,geom,specimen.couplant.pos,specimen.couplant.centerpos)))) |  ...
	  ... % Couplant excitation
	  @(specimen) AddBoundaryCondition(M,specimen,specimen.couplant,[specimen.couplant.tag '_harmonicforceexc'], ...
					   'solidmech_harmonic', ...
					   'impulseforceexcitation', ...
					   @(M,physics,bcobj) ...
					    BuildFaceTotalForceBC(M,geom,physics,specimen.couplant,bcobj, ...
								  specimen.couplant.getfreefaceselection, ...
								  GetOutwardNormal(M,geom,specimen.couplant.pos,specimen.couplant.centerpos)));



% Define a procedure for building the crack. Steps can be sequenced by using
% the pipe (vertical bar | ) character. 
cracktype = 'penny';
bldcrack = @(M,geom,specimen) CreateCrack(M,geom,'crack',specimen, ...
					  { ObtainDCParameter(M,'simulationcrackx','m'), ...
					    ObtainDCParameter(M,'simulationcracky','m'), ...
					    ObtainDCParameter(M,'simulationcrackz','m') }, ...
					  ObtainDCParameter(M,'cracksemimajoraxislen'), ...
					  ObtainDCParameter(M,'cracksemiminoraxislen'), ...
					  [0,1,0], ...
					  [0,0,-1], ...
					  [ .001, .002, .003 ], ...
                                          {'solidmech_multisweep'}, ...
					  dc_dummy_heatingdata_href{1},... % Text file to hold crack heating energies. 
                      cracktype);

% Define a procedure for creating the various physics definitions. Steps can be 
% sequenced by using the pipe (vertical bar | ) character. 
bldphysics = @(M,geom,specimen,flaw) VibroPhysics(M,geom,specimen,flaw,'welderprocess') | @(M) AddXducerContactProbe(M,geom,specimen,couplant_coord(1:3));

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
