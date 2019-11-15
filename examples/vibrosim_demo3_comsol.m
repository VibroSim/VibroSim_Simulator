% function ret = vibrosim_demo3(output_mph_name_href)
%> @brief Here is a third example of how to use BuildVibroModel
%> The commented function declaration is used by processtrak to
%> figure out the parameters to pass. 

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
default_params(M);

% Extract parameters that will be needed below
ObtainDCParameter(M,'staticload_mount');
ObtainDCParameter(M,'xducerforce','N');

% Define a procedure for building the geometry. Steps can be sequenced by using
% the pipe (vertical bar | ) character. 
bldgeom = @(M,geom) CreateRectangularBarSpecimen(M,geom,'specimen') | ...
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
								  {'0','0','-staticload_mount'})) | ...
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
												GetOutwardNormal(M,geom,specimen.couplant.pos,specimen.couplant.centerpos))));


% Define a procedure for building the crack. Steps can be sequenced by using
% the pipe (vertical bar | ) character. 
bldcrack = @(M,geom,specimen) CreateCrack(M,geom,'crack',specimen, ...
					  { ObtainDCParameter(M,'simulationcrackx','m'), ...
					    ObtainDCParameter(M,'simulationcracky','m'), ...
					    ObtainDCParameter(M,'simulationcrackz','m') }, ...
					  ObtainDCParameter(M,'cracksemimajoraxislen'), ...
					  ObtainDCParameter(M,'cracksemiminoraxislen'), ...
					  [0,1,0], ...
					  [0,0,-1], ...
					  [ .001, -30 ; .002, 0 ; .003, 60], ...
                                          {'solidmech_harmonicsweep','solidmech_harmonicburst'}, ...
					  'vibroheatconvert');

% Define a procedure for creating the various physics definitions. Steps can be 
% sequenced by using the pipe (vertical bar | ) character. 
bldphysics = @(M,geom,specimen,flaw) VibroPhysics(M,geom,specimen,flaw,'process');

% Define a procedure to create needed COMSOL result nodes. Generally can pass
% either the wrapped model (M) or unwrapped (model). Pipelining will only work
% with the wrapped model 'M'.
genresults = @(M) VibroResults(M);

% Define a path and filename for saving the generated COMSOL model.
% BuildVibroModel() will create the model's directory if necessary
%savefilename = fullfile(tempdir,sprintf('vibrosim_%s',char(java.lang.System.getProperty('user.name'))),'vibrosim_demo3.mph');

savefilename = output_mph_name_href{1}


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

ret = { { 'generatedmph', { savefilename } } };
