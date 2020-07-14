% function ret = gear_comsol(dc_dest_href,dc_measident_str,dc_dummy_heatingdata_href,dc_spcrayleighdamping_alpha_float,dc_spcrayleighdamping_beta_float,dc_simulationeigsaround_float,dc_excitation_frequency_float,dc_sweep_start_frequency_float,dc_sweep_end_frequency_float,dc_sweep_step_frequency_float,dc_exc_t0_float,dc_exc_t4_float,dc_amplitude_float)
% Here is a third example of how to use BuildVibroModel

[M,model]=InitializeVibroSimScript();

% This code will create a VibroSim model for a gear that is available for
% purchase at McMaster Carr.

% This model is not feature complete.

% PN 5172T45
% 28 Teeth
% Pitch Diameter 3.5 in
% Outer Diameter 3.75 in
% Face Width 1.5 in
% Overall Width 2.38 in
% Hub Diameter 3.06 in
% Hub Width 0.88 in
% Number of teeth 28

meshsizemin = 0.1e-3 ;
meshsize = 3e-3 ; 

AddParamToParamdb(M,'meshsizemin',meshsizemin,'m') ;
AddParamToParamdb(M,'meshsize',meshsize,'m') ;


AddParamToParamdb(M,'simulationneigs',40);  % Number of frequencies for modal analysis to seek out
AddParamToParamdb(M,'simulationeigsaround',dc_simulationeigsaround_float ,'Hz'); % Center frequency for modal analysis

% % simulationtimestart, simulationtimestep, and simulationtimeend specify the time range of the heat flow simulation
AddParamToParamdb(M,'simulationtimestart',dc_exc_t0_float,'s');
AddParamToParamdb(M,'simulationtimestep',0.02,'s');
AddParamToParamdb(M,'simulationtimeend',dc_exc_t4_float,'s'); % .8 seconds after assumed end of vibration

AddParamToParamdb(M,'simulationcameranetd',.022,'K');

% Laser (displacement or velocity detection) coordinates
Hyp = 1.59*25.4e-3; % This is the X location of the vibrometer point on the cross section
Z = 1.5*25.4e-3; % This is the Y location of the vibrometer point on the cross section

% Laser is pointed at area around the hole
laserx = Hyp ;
lasery = 0 ;
laserz = Z ;
AddParamToParamdb(M,'laserx',laserx,'m'); 
AddParamToParamdb(M,'lasery',lasery,'m');
AddParamToParamdb(M,'laserz',laserz,'m');

% Laser (displacement or velocity detection) direction vector
laserdx = 0 ; % cos(toRadians('degrees',0));
laserdy = 0 ; % sin(toRadians('degrees',0));
AddParamToParamdb(M,'laserdx',laserdx);
AddParamToParamdb(M,'laserdy',laserdy);
AddParamToParamdb(M,'laserdz',-1);

% These will be plugged in during processing =======================================
AddParamToParamdb(M,'simulationfreqstart',dc_sweep_start_frequency_float,'Hz');
AddParamToParamdb(M,'simulationfreqstep',dc_sweep_step_frequency_float,'Hz');
AddParamToParamdb(M,'simulationfreqend',dc_sweep_end_frequency_float,'Hz');
AddParamToParamdb(M,'simulationburstfreq',dc_excitation_frequency_float,'Hz');
% ==================================================================================

% NOTE: isolatornlayers removed -- assume something else 
% is doing the multiplication to determine total isolatorthickness
AddParamToParamdb(M,'isolatorthickness',.003,'m');
AddParamToParamdb(M,'isolatorlength',.010,'m');
AddParamToParamdb(M,'isolatorwidth',.0254,'m');

%AddParamToParamdb(M,'isolatormeshtype','TETRAHEDRAL');
AddParamToParamdb(M,'isolatormeshtype','HEXAHEDRAL');
AddParamToParamdb(M,'isolatorfacemethod','FreeQuad');
AddParamToParamdb(M,'isolatormeshsize',.004,'m');
AddParamToParamdb(M,'isolatorsweepelements',6,'');
AddParamToParamdb(M,'isolatorThermalConductivity',.05,'W/m/K');  % this number may not be very meaningful
AddParamToParamdb(M,'isolatorSpecificHeatCapacity',2500,'J/kg/K');  % this number may not be very meaningful

% WARNING: These values do not correspond to any particular material
AddParamToParamdb(M,'isolatormaterial','isocardstock');
AddParamToParamdb(M,'isolatorYoungsModulus',20,'GPa');
AddParamToParamdb(M,'isolatorDensity',800,'kg/m^3');
AddParamToParamdb(M,'isolatorPoissonsRatio',0.3,'');

AddParamToParamdb(M,'couplantthickness',.001,'m');
AddParamToParamdb(M,'couplantlength',.010,'m');
AddParamToParamdb(M,'couplantwidth',.010,'m');

AddParamToParamdb(M,'couplantmeshtype','HEXAHEDRAL');
AddParamToParamdb(M,'couplantfacemethod','FreeQuad');
AddParamToParamdb(M,'couplantmeshsize',.004,'m');
AddParamToParamdb(M,'couplantsweepelements',4,'');

% WARNING: These values do not correspond to any particular material
AddParamToParamdb(M,'couplantmaterial','coupcardstock');
AddParamToParamdb(M,'couplantYoungsModulus',200,'MPa');
AddParamToParamdb(M,'couplantDensity',800,'kg/m^3');
AddParamToParamdb(M,'couplantPoissonsRatio',0.3,'');
AddParamToParamdb(M,'couplantEta',.05,'');
AddParamToParamdb(M,'couplantThermalConductivity',.05,'W/m/K');  % this number may not be very meaningful
AddParamToParamdb(M,'couplantSpecificHeatCapacity',2500,'J/kg/K');  % this number may not be very meaningful
AddParamToParamdb(M,'couplantdashpotcoeff',100e3,'Pa*s/m');

AddParamToParamdb(M,'amplitude',dc_amplitude_float,'V');

AddParamToParamdb(M,'xducercalib',fullfile(fileparts(which('BuildVibroModel')),'..','transducer_data','constant_10micronpervolt_displacementampl.dat'));

CreateTransducerDisplacementVariable(M,fullfile(fileparts(which('BuildVibroModel')),'..','transducer_data','constant_10micronpervolt_displacementampl.dat'),dc_amplitude_float);

CreateCameraNoise(M,'cameranoise',ObtainDCParameter(M,'simulationcameranetd','K'));

% Need the surface emissivity as well
AddParamToParamdb(M,'simulationsurfaceemissivity',1.0);
ObtainDCParameter(M,'simulationsurfaceemissivity');

% Specimen Specific Parameters
% Use the new damping model
AddParamToParamdb(M,'spcmaterialdampingtype','RayleighDamping','');
AddParamToParamdb(M,'spcrayleighdamping_alpha',dc_spcrayleighdamping_alpha_float,'1/s');
AddParamToParamdb(M,'spcrayleighdamping_beta',dc_spcrayleighdamping_beta_float,'s');

% Transducer
couplant_trace_x = 1.25*25.4e-3 ; 
couplant_trace_y = (1.5+0.88)*25.4e-3 ;
couplantx = couplant_trace_x * cos(toRadians('degrees',0));
couplanty = couplant_trace_x * sin(toRadians('degrees',0));
couplantz = couplant_trace_y;

% Welder
% couplant_trace_x = 329.13e-3 ;
% couplant_trace_y =  31.29e-3 ;
% couplantx = couplant_trace_x * cos(toRadians('degrees',10*7.2+90));
% couplanty = couplant_trace_x * sin(toRadians('degrees',10*7.2+90));
% couplantz = couplant_trace_y;

isolator_coords=[];
couplant_coord = [ couplantx couplanty couplantz NaN]; 

% crack_coord=[36.621e-3  109.835e-3  -.0078091]; 
Hyp = 1.605*25.4e-3; % This is the X location of the vibrometer point on the cross section
Z = 1.5*25.4e-3; % This is the Y location of the vibrometer point on the cross section
crackx = Hyp - 0.1e-3
cracky = 0
crackz = Z - 0.1e-3
crack_coord=[ crackx cracky crackz ];

% Crack semiminor axis (through the depth of the material)
crack_minor_axis= [ 0 0 -1 ] ; % crack_coord./norm(crack_coord)-toRadians('degrees',30); % This is for a crack that is 30 degrees offset from the radial

% Find the semimajor axis
major_axis_p1 = [crackx cracky crackz]; % Start with the vector point to the crack center
major_axis_p2 = [ -1 0 0 ] ;

% crack_minor_axis = (csm_axis_p2-csm_axis_p1)./norm(csm_axis_p2-csm_axis_p1);
crack_major_axis = (major_axis_p2-major_axis_p1)./norm(major_axis_p2-major_axis_p1);

crack_semimajor_len=2e-3;
crack_semiminor_len=1e-3;

crack_closure_steps=6;  % For some reason we run into geometry failures (negative domain volume!)
                        % If we give this too many steps (?)
crack_closure_step=crack_semimajor_len/crack_closure_steps;
crack_closure=[ (crack_closure_step:crack_closure_step:crack_semimajor_len).' ];

% Crack type --> 'penny' or 'through'
cracktype = 'penny';

gear_cross_section_file = fullfile(pwd,'vibrosim_demo5_geometry','gear_cross_section.dxf');
bldgeom = @(M,geom) gear_geometry(M,geom,'specimen',gear_cross_section_file,meshsizemin,meshsize)  | ...
        @(specimen) AttachThinCouplantIsolators(M,geom,specimen, ...
					couplant_coord, ...
					isolator_coords) | ...
        ... % Force condition on couplant
           @(specimen) AddBoundaryCondition(M,specimen,specimen.couplant,[specimen.couplant.tag '_xducerforce'], ...
                 'solidmech_static', ...
                 'staticloading', ...
                 @(M,physics,bcobj) ...
                  BuildFaceTotalForceBC(M,geom,physics,specimen.couplant,bcobj, ...
                            specimen.couplant.getfreefaceselection, ...
                            ... % xducerforce is multiplied by outward normal 
                                        ... % on the couplant, which is equivalent to        the
                                        ... % inward normal on the specimen
                            MultiplyScalarStrByNumericVec('xducerforce', ...
                                          GetOutwardNormal(M,geom,specimen.couplant.pos,specimen.couplant.centerpos))));
         

% bldcrack = [];
bldcrack = @(M,geom,specimen) CreateCrack(M,geom,'crack',specimen, ...
                    crack_coord, ...
                    crack_semimajor_len, ...
                    crack_semiminor_len, ...
                    crack_major_axis, ...
                    crack_minor_axis, ...
                    crack_closure, ...
                    {'solidmech_harmonicsweep', 'solidmech_harmonicburst'},...
					dc_dummy_heatingdata_href{1}, ...; % Text file to hold crack heating energies. 
                    cracktype); 

bldphysics = @(M,geom,specimen,flaw) VibroPhysics(M,geom,specimen,flaw,'broadbandprocess');
% bldphysics = [] ;

genresults = @(M) VibroResults(M);
% genresults = [] ;

% savefilename = fullfile('/','tmp',sprintf('AnalyzeF101_%s.mph',char(java.lang.System.getProperty('user.name'))));
savefilename = fullfile(dc_dest_href{1},[ dc_measident_str '.mph']);

% Given the model wrapper, the procedures for building the geometry, the crack, the physics,
% and for storing the results, and a file name for saving, have COMSOL build the model.x
BuildVibroModel(M,...
		bldgeom, ...
		bldcrack, ...
		bldphysics, ...
		genresults, ...
		savefilename);

ret = { { 'dc:model_comsol', { savefilename } } };
