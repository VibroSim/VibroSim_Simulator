% function ret = cantilever_example_viscousdamping_comsol(dc_dest_href,dc_measident_str,dc_dummy_heatingdata_href,dc_xducerforce_float,dc_spcmaterial_str,dc_YoungsModulus_float, dc_YieldStrength_float, dc_PoissonsRatio_float, dc_Density_float,dc_spcviscousdamping_float,dc_mountdamping_scalefactor_float,dc_mountstiffness_scalefactor_float,dc_baseline_mountstiffness_float,dc_limiting_mountdamping_float,dc_couplantx_float, dc_couplanty_float, dc_couplantz_float, dc_couplantangle_float, dc_exc_t0_float, dc_exc_t4_float)
%> @brief Here is a variant on the cantilever example using viscous
%> material damping and a model of radiative damping at the cantilever mount
%>
%> Note that because the radiative damping is frequency-dependent the 
%> modal analysis may not capture it properly except right near the 
%> linearization point of the modal solver (settable with 
%> simulationeigsaround parameter. As a result the synthetic spectrum
%> may not be much good and manual tweaking of the freqbands may 
%> be required. This would usually be done by running the multisweep
%> analysis with much larger step sizes, so as to use the mode widths
%> to determine q-factors and resonant decay times. 
%>
%> for this cantilever_example_viscousdamping one set of verified 
%> segment parameters are believed to be (CHECK .mph file from run!): 
%>seg1 start 0.000
%>seg1 step 0.476837158203125
%>seg1 end 74.3865966796875
%>seg2 start 73.909759521484375
%>seg2 step 7.62939453125 
%>seg2 end 2660.27450561523438
%>seg3 start 2662.65869140625
%>seg3 step 61.03515625
%>seg3 end 23231.50634765625
%>seg4 start 23236.7515563964844
%>seg4 step 488.28125
%>seg4 end 999999.523162841797

% NOTE: The first line (commented out function definition)
%       is important because it tells processtrak which parameters
%       to pass
%
% NOTE: The name of this script should end with "_comsol.m" as that
%       makes processtrak run it with COMSOL MATLAB rather than
%       bare MATLAB>


%
% NOTE: Make your window as wide as possible when viewing this file!



% InitializeVibroSimScript() connects to COMSOL and initializes and 
% returns the wrapped model variable (M) and the unwrapped node (model). 
[M,model]=InitializeVibroSimScript();

% Add additional parameters...

% spcviscousdamping is used by cantilever_example_viscousdamping_comsol.. but not 
% cantilever_example_comsol.m because only the former sets spcmaterialdampingtype to ViscousSamping

%AddParamToParamdb(M,'spcviscousdamping',2e3,'N*s'); % Value used after 11/06/19
AddParamToParamdb(M,'spcviscousdamping',dc_spcviscousdamping_float,'N*s'); 

AddParamToParamdb(M,'mountdamping_scalefactor',dc_mountdamping_scalefactor_float); % Factor by which the mount damping is larger than calculated from theory ... 2.5e-5 for 11/7/19
AddParamToParamdb(M,'mountstiffness_scalefactor',dc_mountstiffness_scalefactor_float); % Factor by which the mount stiffness is larger than calculated from theory... 2.5e-5 for 11/7/19

% Baseline mountstiffness and limiting_mountdamping are to help make the frequency-dependent mount stiffness and mount 
% damping closer to linear to help with the eigenvalue (modal) analysis. 
% Generally select them by adding frequency domain 1D plots of the cantilever_mount_dampingperarea and cantilever_mount_stiffnessperarea
% and making them low and high enough respectively not to mess too much with the lowest order resonance of interest. 
AddParamToParamdb(M,'baseline_mountstiffness',dc_baseline_mountstiffness_float,'N/m^3');
AddParamToParamdb(M,'limiting_mountdamping',dc_limiting_mountdamping_float,'N*s/m^3');



 
% Call a function that sets various parameters to be used by the model. 
cantilever_modeling_params(M);


% Set COMSOL parameters from experiment log, passed as parameters
AddParamToParamdb(M,'xducerforce',dc_xducerforce_float,'N');
AddParamToParamdb(M,'spcmaterial',dc_spcmaterial_str);
AddParamToParamdb(M,'spcYoungsModulus',dc_YoungsModulus_float,'Pa');
AddParamToParamdb(M,'spcPoissonsRatio',dc_PoissonsRatio_float,'');
AddParamToParamdb(M,'spcDensity',dc_Density_float,'kg/m^3');

% simulationtimestart, simulationtimestep, and simulationtimeend specify the time range of the heat flow simulation
AddParamToParamdb(M,'simulationtimestart',dc_exc_t0_float,'s');
AddParamToParamdb(M,'simulationtimestep',0.02,'s');
AddParamToParamdb(M,'simulationtimeend',dc_exc_t4_float+0.8,'s'); % .8 seconds after assumed end of vibration

AddParamToParamdb(M,'spcmaterialdampingtype','RayleighDamping');
AddParamToParamdb(M,'spcrayleighdamping_alpha',dc_spcrayleighdamping_alpha_float,'s^-1');
AddParamToParamdb(M,'spcrayleighdamping_beta',dc_spcrayleighdamping_beta_float,'s');

%                 x          y         z      angle
%couplant_coord=[ .245,       .025,        0,      NaN    ];
couplant_coord=[ dc_couplantx_float,  dc_couplanty_float,  dc_couplantz_float,  dc_couplantangle_float   ];
isolator_coords=[];
% ( NaN for angle causes it to create a circular isolator)

  

% Extract needed parameters...

ObtainDCParameter(M,'xducerforce','N');


% Fixed End springiness and damping characteristics 
% approximated based on:
% Gazetas, G., and Dobry, R. Simple Radiation Damping Model for Piles and Footings, J. Eng. Mech 1984 110(6) 937-956

% According to the paper, 
% c = rho*v*A*real(-i*besselh(1,2,a)./besselh(0,2,a)) where rho is density, v is substrate wavespeed, and A is area.
% where a = 2*pi*f*B/V where A=wB   (Eq. 16,17)
% this function gets moderately large for very small a and asymptotes to 1.0 for large a... a represents f*width/wavespeed.... so at very low frequencies it can get very big....
% Alternatively, the Lysimer model (Eq. 1) in the paper is frequency independent... cv = rho*A*3.4*Vs/(pi*(1-nu))

% For this 2 inch wide and 1/4 inch thick cantilever and low order bending modes @ ~ 100 Hz: 
%   * In this mode based on the diagrams in the paper the relevant B is half the thickness, giving A=2B=thickness=6.35 mm
%      with rho = 8050 kg/m^3 and V = 3100 m/s for a vertical shear wave in the steel substrate, nu~0.3 and V = 5790 m/s for longitudinal wave
%      and Vla = 4790 m/s
%   * Lysimer model: cv = 12400 N/(m/s)
%   * Eq. 16 then gives, for f=20 kHz
%      Eq. 16: a = 0.129,  the real part of the Henkel function quotient is 1.72 and c = (273000 N/(m/s))/m. For the 
%      2 inch wide cantilever this corresponds to a total damping coefficient of 13800 N/(m/s)... This is for shear waves
%      propagating 'vertically', meaning normal to the substrate surface. 
%   * Alternatively if the substrate is a plate we would use the pile model generating Sh waves and p-waves
%      Adjusting eq. 20 for our rectangular support in that lowest resonance, 
%      we get ch = 2*(thickness)*(substrate_thickness)*rho*Vs*(real part of henkel quotient with Vs) +
%                  + 2*(width)*substrate_thickness)*rho*Vs*(3.4/(pi*(1-nu)))*(real part of henkel quotient with Vla)
%          The real part of henkel quotient with Vs would come from omega*(6.35mm/2)/Vs or .129 
%          The real part of henkel quotient with Vla would come from omega*(50.8mm/2)/Vla or .666
%      for 1/2" thick steel substrate this gives
%      ch = 4025 N/(m/s) * 1.72        + 32200 N/(m/s) * 1.13 ... or
%      ch = 43300 N/(m/s)  -- damping coefficient with respect to horizontal motions
%      ... this is proportional to depth which was 1/2". Expect effective depth to be related to thickness of cantilever (1/4"), so maybe expect half this
%      ... If we assume damping comes from all of these sources, 
%      get c = 13800 + 43300 = 60100 N/(m/s)
%   * Should also give end some springiness for DC stability
%     Spring in parallel with dashpot = (jw/k) || (1/c)  -> (jw/ck)/(jw/k + 1/c) -> (jw)/(jwc +k)
%   want k << jwc at resonance.... assuming w = 2*pi*100 Hz and c of 60000 N/(m/s), k << 37.7e6 N/m 
% ... use k = 1 N/um or 1e6 N/m

% plot(a,real(-i*besselh(1,2,a)./besselh(0,2,a)))

% ***!!!NOTE: Needed to adjust modal solver to get resonances beyond 20 kHz!!!***

% Frequency-dependent parameterized calculation of footing model from above
spclength = ObtainDCParameter(M,'spclength','m');
spcwidth = ObtainDCParameter(M,'spcwidth','m');  % 50.8 mm for cantilever
spcthickness = ObtainDCParameter(M,'spcthickness','m');  % 6.35 mm for cantilever

mountdamping_scalefactor = ObtainDCParameter(M,'mountdamping_scalefactor');
mountstiffness_scalefactor = ObtainDCParameter(M,'mountstiffness_scalefactor');

% The frequency-dependent damping defined here is highly nonlinear near the origin (zero frequency)
% The eigenvalue solver used for modal analysis cannot handle nonlinear frequency dependence  -- it linearizes. 
% So we do two things to correct and accommodate the frequency-dependent damping: 
%   1. The mount stiffness approaches zero as freq->0, We add in a baseline stiffness to prevent it from getting too small. 
%   2. The mount damping gets large as freq->0. We parallel in a limiting damping coefficient to prevent it from getting too large. 
%   3. We ask the eigenvalue solver to do its linearization about a relevant frequency (i.e. not zero).  CreateVibroModal() uses the 'simulationeigsaround' parameter for this. 
baseline_mountstiffness = ObtainDCParameter(M,'baseline_mountstiffness');
limiting_mountdamping = ObtainDCParameter(M,'limiting_mountdamping');

%substrate_V = 5790; % Pressure velocity for steel substrate, m/s
substrate_Vs=3100;  % Shear velocity for steel substrate, m/s
substrate_nu=0.3; % poisson's ratio for steel substrate
substrate_Vla=3.4*substrate_Vs/(pi*(1-substrate_nu)); % Vla (Lysmer's wave velocity)
substrate_rho = 8050; % density of steel substrate, kg/m^3

% Hankel functions of the second kind as COMSOL expressions
HankelSecondKind = @ (alpha,x) [ '(besselj(' alpha ',' x ')' ' - i*bessely(' alpha ',' x '))'];

% real parts of Hankel quotients from the paper
hankelquotient_Vs = [ '(-i*(' HankelSecondKind('1',[ '(pi/4)*2*pi*(abs(real(freq))+.001[Hz])*((' spcthickness ')/2.0)/(' num2str(substrate_Vs) '[m/s]' ')']) '/' HankelSecondKind('0',[ '(pi/4)*2*pi*(abs(real(freq))+.001[Hz])*((' spcthickness ')/2.0)/(' num2str(substrate_Vs) '[m/s]' ')']) '))' ];
hankelquotient_Vla = [ '(-i*(' HankelSecondKind('1',[ '(pi/4)*2*pi*(abs(real(freq))+.001[Hz])*((' spcwidth ')/2.0)/(' num2str(substrate_Vla) '[m/s]' ')']) '/' HankelSecondKind('0',[ '(pi/4)*2*pi*(abs(real(freq))+.001[Hz])*((' spcwidth ')/2.0)/(' num2str(substrate_Vla) '[m/s]' ')']) '))'];

% Calculation of ch
%      we get ch = 2*(thickness)*(substrate_thickness)*rho*Vs*(real part of henkel quotient with Vs) +
%                  + 2*(width)*substrate_thickness)*rho*Vs*(3.4/(pi*(1-nu)))*(real part of henkel quotient with Vla) (the Vs*3.4/(pi*(1-nu)) is equiv to Vla)

% Treat effective substrate thickness as twice specimen thickness
ch = [ 'abs(2.0*(' spcthickness ')*2.0*(' spcthickness ')*' num2str(substrate_rho) '[kg/m^3]' '*' num2str(substrate_Vs) '[m/s]' '* real(' hankelquotient_Vs ')+' ...
       '2.0*(' spcwidth ')*2.0*(' spcthickness ')*' num2str(substrate_rho) '[kg/m^3]' '*' num2str(substrate_Vla) '[m/s]' '* real(' hankelquotient_Vla '))'];

a = ['(2*pi*(abs(real(freq))+.001[Hz])*((' spcthickness ')/2)/' num2str(substrate_Vs) '[m/s]' ')'];
c_eq16 = [ 'abs(' num2str(substrate_rho) '[kg/m^3]' '*' num2str(substrate_Vs) '[m/s]' '*(' spcthickness ')*(' spcwidth ')*real(-i*(' HankelSecondKind('1',a) ')/(' HankelSecondKind('0',a) ')))' ];

ctotal_perarea = ['((((' c_eq16 ') + (' ch '))*(' mountdamping_scalefactor '))/((' spcthickness ')*(' spcwidth ')))']; % kg/(m^2*s)  ... equiv to N*s/m^3 = (kg*m/s^2)*s/m^3 = kg/(s*m^2)  CHECK 

% implement limiting_mountdamping (see above): parallel combination of ctotal_perarea with limiting_mountdamping constant
climited = ['(1/(' '1/(' ctotal_perarea ')' '+' '1/(' limiting_mountdamping ')'  '))'];

CreateVariable(M,'cantilever_mount_dampingperarea',climited);

% Gazetas et al. also gives us a way to predict mounting stiffnesses. 
% Specifically, the c's they give are real(F0/udot0) = real(F0/iwu0)
% ... The stiffness K = real(F0/u0) = real(F0*iw/iwu0) = -imag(F0*w/iwu0)
%  =-w*imag(F0/iwu0)
% So we do the same calculation as for damping coefficient but
% take imaginary part instead of real and multiply by -w = -2*pi*f

kh = [ 'abs((-2*pi*real(freq))*(2.0*(' spcthickness ')*2.0*(' spcthickness ')*' num2str(substrate_rho) '[kg/m^3]' '*' num2str(substrate_Vs) '[m/s]' '* imag(' hankelquotient_Vs ')+' ...
       '2.0*(' spcwidth ')*2.0*(' spcthickness ')*' num2str(substrate_rho) '[kg/m^3]' '*' num2str(substrate_Vla) '[m/s]' '* imag(' hankelquotient_Vla ')))'];


k_eq16 = [ 'abs(' '-2*pi*real(freq)*' num2str(substrate_rho) '[kg/m^3]' '*' num2str(substrate_Vs) '[m/s]' '*(' spcthickness ')*(' spcwidth ')*imag(-i*(' HankelSecondKind('1',a) ')/(' HankelSecondKind('0',a) ')))' ];

ktotal_perarea = ['((' baseline_mountstiffness ') + (((' k_eq16 ') + (' kh '))*(' mountstiffness_scalefactor '))/((' spcthickness ')*(' spcwidth ')))']; 

CreateVariable(M,'cantilever_mount_stiffnessperarea',ktotal_perarea);

% Damping tensor for 'fixed' and of cantilever; represents absorptivity of mount 
% units: (N*s/m)/m^2   ... should be a 9 element cell array representing elements of 3x3 tensor
%FixedEnd_DampPerArea = to_cellstr_array([ 60100/(6.35e-3*50.8e-3), 0, 0,
%					  0, 60100/(6.35e-3*50.8e-3), 0,
%					  0, 0, 60100/(6.35e-3*50.8e-3)]);
FixedEnd_DampPerArea = { 'cantilever_mount_dampingperarea', '0', '0', ...
			 '0', 'cantilever_mount_dampingperarea', '0', ...
			 '0', '0', 'cantilever_mount_dampingperarea' };


% Stiffness tensor for 'fixed' and of cantilever; represents springiness of mount 
% units: (N/m)/m^2   ... should be a 9 element cellstr array representing elements of 3x3 tensor
%% ... XXX Just make stiffness & damping the same in y and z axes for now... for x axis use stiffness of 1cm steel over a fixed BC to support boundary against bending moments
%FixedEnd_k_A  = to_cellstr_array([ 200e9/1e-2, 0, 0, 
%				   0, 200e9/1e-2, 0, 
%				   0, 0, 200e9/1e-2]);

%				   0, 1e6/(6.35e-3*50.8e-3), 0, 
%				   0, 0, 1e6/(6.35e-3*50.8e-3)]);
FixedEnd_k_A = { 'cantilever_mount_stiffnessperarea', '0', '0', ...
		 '0', 'cantilever_mount_stiffnessperarea', '0', ...
		 '0', '0', 'cantilever_mount_stiffnessperarea' };



% We just make the stiffness complex and forget about the Damping entry in COMSOL

% For spring x=F/k_A
% for damping jwx = F/c  or x=F/jwc
% We have x = F/k_A + F/jwc = F*(1/k_A  + 1/jwc)
% k_A_eff = 1/(1/k_A _ 1/jwc)

FixedEnd_k_A_combined = { ['(' FixedEnd_k_A{1} ') + (i*2*pi*freq*cantilever_mount_dampingperarea)'], '0', '0', ...
			  '0', ['(' FixedEnd_k_A{5} ') + (i*2*pi*freq*cantilever_mount_dampingperarea)'], '0', ...
			  '0', '0', ['(' FixedEnd_k_A{9} ') + (i*2*pi*freq*cantilever_mount_dampingperarea)']};

% !!!*** NOTE: This damping is only compatible with frequency domain...




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
					    ... %BuildFaceFixedBC(M,geom,physics,specimen,bcobj,...
					    ... %                  @(M,geom,spec) GetBlockFace(M,geom,specimen,[-1,0,0]))) | ...
					    ... BuildFaceSpringFoundationBC(M,geom,physics,specimen,bcobj,...
					    				    @(M,geom,spec) GetBlockFace(M,geom,specimen,[-1,0,0]), ...
					    				    ... %FixedEnd_k_A, FixedEnd_DampPerArea)) |...
					     				    FixedEnd_k_A_combined)) |...
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
