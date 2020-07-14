function specimen=gear_geometry(M,geom,tag,crosssect_dxf_file,meshsizemin,meshsize)

specimen=ModelWrapper(M,tag);

% Create work plane for cross section. Import DXF file
CreateWrappedProperty(M,specimen,'crosssect',[tag '_crosssect'],geom.node,'WorkPlane');
specimen.crosssect.node.set('unite','on');
specimen.crosssect.node.set('quickplane','yz');
specimen.crosssect.node.set('createselection','on');

% Create import node
CreateWrappedProperty(M,specimen,'crosssect_import',[tag '_crosssect_import'],specimen.crosssect.node.geom,'Import');
specimen.crosssect_import.node.set('type','dxf');
specimen.crosssect_import.node.set('filename',crosssect_dxf_file);
specimen.crosssect_import.node.importData;

% Create scaling node to convert dxf file (in mm) to our units (m)
CreateWrappedProperty(M,specimen,'crosssect_scale',[tag '_crosssect_scale'],specimen.crosssect.node.geom,'Scale');
specimen.crosssect_scale.node.selection('input').set({specimen.crosssect_import.tag});
specimen.crosssect_scale.node.set('isotropic','1e-3');

% Revolve the part
CreateWrappedProperty(M,specimen,'revolve',[tag '_revolve'],geom.node,'Revolve');
specimen.revolve.node.set('angle2','90');
specimen.revolve.node.selection('input').set({specimen.crosssect.tag});
specimen.revolve.node.set('createselection','on');


%% Cutting the teeth
% Create work plane for gear teeth. Import DXF file
CreateWrappedProperty(M,specimen,'teeth',[tag '_teeth'],geom.node,'WorkPlane');
specimen.teeth.node.set('unite','on');
% specimen.teeth.node.set('quickplane','xy');
specimen.teeth.node.set('quickz','-1e-3'); 
specimen.teeth.node.set('createselection','on');

% CreateWrappedProperty(M,specimen,'holeplane',[tag '_holeplane'],geom.node,'WorkPlane');

% Create import node
gear_tooth_file = fullfile(pwd,'vibrosim_demo5_geometry','gear_tooth.dxf');
CreateWrappedProperty(M,specimen,'teeth_import',[tag '_teeth_import'],specimen.teeth.node.geom,'Import');
specimen.teeth_import.node.set('type','dxf');
specimen.teeth_import.node.set('filename',gear_tooth_file);
specimen.teeth_import.node.importData;

% Create scaling node to convert dxf file (in mm) to our units (m)
CreateWrappedProperty(M,specimen,'teeth_scale',[tag '_teeth_scale'],specimen.teeth.node.geom,'Scale');
specimen.teeth_scale.node.selection('input').set({specimen.teeth_import.tag});
specimen.teeth_scale.node.set('isotropic','1e-3');

% Array of thetas
% There are 28 teeth total, 12.85714 degree spacing
% A quarter of a gear has 8 teeth cut into it.
theta = [0:7]*12.85714;
H = 1.605*25.4*1e-3;

% Create a move node and rotate it [1]
CreateWrappedProperty(M,specimen,'tooth_move1',[tag '_tooth_move1'],specimen.teeth.node.geom,'Move');
specimen.tooth_move1.node.selection('input').set({specimen.teeth_scale.tag});
specimen.tooth_move1.node.set('keep',true);
specimen.tooth_move1.node.set('disply',num2str(H)) ;
CreateWrappedProperty(M,specimen,'tooth_rotate1',[tag '_tooth_rotate1'],specimen.teeth.node.geom,'Rotate');
specimen.tooth_rotate1.node.selection('input').set({specimen.tooth_move1.tag});
specimen.tooth_rotate1.node.set('rot',theta(1));

% Create a move node and rotate it [2]
CreateWrappedProperty(M,specimen,'tooth_move2',[tag '_tooth_move2'],specimen.teeth.node.geom,'Move');
specimen.tooth_move2.node.selection('input').set({specimen.teeth_scale.tag});
specimen.tooth_move2.node.set('keep',true) ;
specimen.tooth_move2.node.set('disply',num2str(H)) ;
CreateWrappedProperty(M,specimen,'tooth_rotate2',[tag '_tooth_rotate2'],specimen.teeth.node.geom,'Rotate');
specimen.tooth_rotate2.node.selection('input').set({specimen.tooth_move2.tag});
specimen.tooth_rotate2.node.set('rot',theta(2));

% Create a move node and rotate it [3]
CreateWrappedProperty(M,specimen,'tooth_move3',[tag '_tooth_move3'],specimen.teeth.node.geom,'Move');
specimen.tooth_move3.node.selection('input').set({specimen.teeth_scale.tag});
specimen.tooth_move3.node.set('keep',true) ;
specimen.tooth_move3.node.set('disply',num2str(H)) ;
CreateWrappedProperty(M,specimen,'tooth_rotate3',[tag '_tooth_rotate3'],specimen.teeth.node.geom,'Rotate');
specimen.tooth_rotate3.node.selection('input').set({specimen.tooth_move3.tag});
specimen.tooth_rotate3.node.set('rot',theta(3));

% Create a move node and rotate it [4]
CreateWrappedProperty(M,specimen,'tooth_move4',[tag '_tooth_move4'],specimen.teeth.node.geom,'Move');
specimen.tooth_move4.node.selection('input').set({specimen.teeth_scale.tag});
specimen.tooth_move4.node.set('keep',true) ;
specimen.tooth_move4.node.set('disply',num2str(H)) ;
CreateWrappedProperty(M,specimen,'tooth_rotate4',[tag '_tooth_rotate4'],specimen.teeth.node.geom,'Rotate');
specimen.tooth_rotate4.node.selection('input').set({specimen.tooth_move4.tag});
specimen.tooth_rotate4.node.set('rot',theta(4));

% Create a move node and rotate it [5]
CreateWrappedProperty(M,specimen,'tooth_move5',[tag '_tooth_move5'],specimen.teeth.node.geom,'Move');
specimen.tooth_move5.node.selection('input').set({specimen.teeth_scale.tag});
specimen.tooth_move5.node.set('keep',true) ;
specimen.tooth_move5.node.set('disply',num2str(H)) ;
CreateWrappedProperty(M,specimen,'tooth_rotate5',[tag '_tooth_rotate5'],specimen.teeth.node.geom,'Rotate');
specimen.tooth_rotate5.node.selection('input').set({specimen.tooth_move5.tag});
specimen.tooth_rotate5.node.set('rot',theta(5));

% Create a move node and rotate it [6]
CreateWrappedProperty(M,specimen,'tooth_move6',[tag '_tooth_move6'],specimen.teeth.node.geom,'Move');
specimen.tooth_move6.node.selection('input').set({specimen.teeth_scale.tag});
specimen.tooth_move6.node.set('keep',true) ;
specimen.tooth_move6.node.set('disply',num2str(H)) ;
CreateWrappedProperty(M,specimen,'tooth_rotate6',[tag '_tooth_rotate6'],specimen.teeth.node.geom,'Rotate');
specimen.tooth_rotate6.node.selection('input').set({specimen.tooth_move6.tag});
specimen.tooth_rotate6.node.set('rot',theta(6));

% Create a move node and rotate it [7]
CreateWrappedProperty(M,specimen,'tooth_move7',[tag '_tooth_move7'],specimen.teeth.node.geom,'Move');
specimen.tooth_move7.node.selection('input').set({specimen.teeth_scale.tag});
specimen.tooth_move7.node.set('keep',true) ;
specimen.tooth_move7.node.set('disply',num2str(H)) ;
CreateWrappedProperty(M,specimen,'tooth_rotate7',[tag '_tooth_rotate7'],specimen.teeth.node.geom,'Rotate');
specimen.tooth_rotate7.node.selection('input').set({specimen.tooth_move7.tag});
specimen.tooth_rotate7.node.set('rot',theta(7));

% Create a move node and rotate it [8]
CreateWrappedProperty(M,specimen,'tooth_move8',[tag '_tooth_move8'],specimen.teeth.node.geom,'Move');
specimen.tooth_move8.node.selection('input').set({specimen.teeth_scale.tag});
specimen.tooth_move8.node.set('disply',num2str(H)) ;
CreateWrappedProperty(M,specimen,'tooth_rotate8',[tag '_tooth_rotate8'],specimen.teeth.node.geom,'Rotate');
specimen.tooth_rotate8.node.selection('input').set({specimen.tooth_move8.tag});
specimen.tooth_rotate8.node.set('rot',theta(8));

% Union the cutting tool
CreateWrappedProperty(M,specimen,'teeth_union',[tag '_teeth_union'],specimen.teeth.node.geom,'Union');
specimen.teeth_union.node.selection('input').set({specimen.tooth_rotate1.tag, ...
                                                  specimen.tooth_rotate2.tag, ...
                                                  specimen.tooth_rotate3.tag, ...
                                                  specimen.tooth_rotate4.tag, ...
                                                  specimen.tooth_rotate5.tag, ...
                                                  specimen.tooth_rotate6.tag, ...
                                                  specimen.tooth_rotate7.tag, ...
                                                  specimen.tooth_rotate8.tag});



CreateWrappedProperty(M,specimen,'teeth_cutting_tool',[tag '_teeth_cutting_tool'],geom.node,'Extrude');
specimen.teeth_cutting_tool.node.set('unite',false);
specimen.teeth_cutting_tool.node.setIndex('distance',to_string(0.1),0);
specimen.teeth_cutting_tool.node.set('crossfaces',false);
specimen.teeth_cutting_tool.node.set('face','all');
specimen.teeth_cutting_tool.node.selection('input').set({specimen.teeth.tag});
specimen.teeth_cutting_tool.node.set('createselection','on');
specimen.teeth_cutting_tool.node.set('inputhandling','remove');


% % subtract teeth tool extrusion
CreateWrappedProperty(M,specimen,'quartergear',[tag '_quartergear'],geom.node,'Difference');
specimen.quartergear.node.selection('input').set({specimen.revolve.tag});
specimen.quartergear.node.selection('input2').set({specimen.teeth_cutting_tool.tag});
specimen.quartergear.node.set('createselection','on');

% mirror the part about the xz plane, then about the yz plane
CreateWrappedProperty(M,specimen,'mirrorxz',[tag '_mirrorxz'],geom.node,'Mirror');
specimen.mirrorxz.node.set('keep',true);
specimen.mirrorxz.node.set('axis', [0,-1,0]);
specimen.mirrorxz.node.set('pos', [1,0,1]);
specimen.mirrorxz.node.selection('input').set({specimen.quartergear.tag});
specimen.mirrorxz.node.set('createselection','on');

CreateWrappedProperty(M,specimen,'mirroryz',[tag '_mirroryz'],geom.node,'Mirror');
specimen.mirroryz.node.set('keep',true);
specimen.mirroryz.node.set('axis', [-1,0,0]);
specimen.mirroryz.node.selection('input').set({specimen.quartergear.tag,specimen.mirrorxz.tag});
specimen.mirroryz.node.set('createselection','on');




%%% IF SKIPPING THE BLADE VANES, UNCOMMENT THE FOLLOWING LINES
% Form a union and place it as the node of specimen
BuildWrappedModel(M,specimen,geom.node,'Union');
specimen.node.selection('input').set({specimen.quartergear.tag,specimen.mirrorxz.tag,specimen.mirroryz.tag});
specimen.node.set('createselection','on');

% Add property to represent specimen or union of specimen with anything else 
% -- such as crack, thin contactors, etc. -- that needs to be added later
addprop(specimen,'specimenunion');
specimen.specimenunion=specimen;

% Set repair tolerance -- gear is complicated enough, esp. with 
% couplant/contactors that we get errors otherwise
geom.node.feature('fin').set('repairtol', '1.0e-9');

geom.node.run;

% Create material
CreateWrappedProperty(M,specimen,'gearmaterial','gearmaterial',M.node.material);
specimen.gearmaterial.node.label(sprintf('1020 Steel (%s)',specimen.gearmaterial.tag));

SetMaterialProperty(M,specimen.gearmaterial,'heatcapacity','486[J/(kg*K)]');
SetMaterialProperty(M,specimen.gearmaterial,'thermalconductivity','51.9[W/(m*K)]');
SetMaterialProperty(M,specimen.gearmaterial,'youngsmodulus', '205e9[Pa]');
SetMaterialProperty(M,specimen.gearmaterial,'poissonsratio', '0.29');
SetMaterialProperty(M,specimen.gearmaterial,'density', '7870[kg/m^3]');

% Apply material 
specimen.applymaterial=BuildLater(M,[specimen.tag '_applymaterial'],...
				    'applymaterial', ...
				    @(M,obj) ...
				    ReferenceNamedMaterial(M,geom,specimen,specimen.gearmaterial.tag,@(M,geom,obj) GetAutomaticSelectionEntities(M,geom,specimen,'dom')));
				

specimen.mesh=BuildLater(M,[specimen.tag '_mesh'],...
			 'meshbuilder',...
			 @(M,mesh,obj) MeshGear(M,geom,mesh,specimen,obj,meshsizemin,meshsize));
