function[specimen]=gear_geometry(M,geom,tag,meshsizemin,meshsize)

%Create a bear model wrapper with no content at the beginning.
specimen=ModelWrapper(M,tag);

% Create import node
CreateWrappedProperty(M,specimen,'imported_geometry',[tag '_imported_geometry'],geom.node,'Import');
specimen.imported_geometry.node.set("type","cad");
specimen.imported_geometry.node.set("filename", "/tmp/al_right_angle_origin.STEP");
specimen.imported_geometry.node.importData

% Form a union and place it as the node of specimen
BuildWrappedModel(M,specimen,geom.node,'Union');
specimen.node.selection('input').set({specimen.imported_geometry.tag});
specimen.node.set('createselection','on');

% Add property to represent specimen or union of specimen with anything else 
% -- such as crack, thin contactors, etc. -- that needs to be added later
addprop(specimen,'specimenunion');
specimen.specimenunion=specimen;

% Set repair tolerance -- gear is complicated enough, esp. with 
% couplant/contactors that we get errors otherwise
geom.node.feature('fin').set('repairtol', '1.0e-9');

geom.node.run;

% Create material for the specimen that was imported.
CreateWrappedProperty(M,specimen,'importmaterial','importmaterial',M.node.material);
specimen.importmaterial.node.label(sprintf('6061 Aluminum (%s)',specimen.importmaterial.tag));

SetMaterialProperty(M,specimen.importmaterial,'heatcapacity','896[J/(kg*K)]');
SetMaterialProperty(M,specimen.importmaterial,'thermalconductivity','167[W/(m*K)]');
SetMaterialProperty(M,specimen.importmaterial,'youngsmodulus', '68.9e9[Pa]');
SetMaterialProperty(M,specimen.importmaterial,'poissonsratio', '0.33');
SetMaterialProperty(M,specimen.importmaterial,'density', '2700[kg/m^3]');
%from: http://www.matweb.com/search/datasheet_print.aspx?matguid=1b8c06d0ca7c456694c7777d9e10be5b

% Apply material 
specimen.applymaterial=BuildLater(M,[specimen.tag '_applymaterial'],...
				    'applymaterial', ...
				    @(M,obj) ...
				    ReferenceNamedMaterial(M,geom,specimen,specimen.importmaterial.tag,@(M,geom,obj) GetAutomaticSelectionEntities(M,geom,specimen,'dom')));
				

specimen.mesh=BuildLater(M,[specimen.tag '_mesh'],...
			 'meshbuilder',...
			 @(M,mesh,obj) MeshAngle(M,geom,mesh,specimen,obj,meshsizemin,meshsize));
