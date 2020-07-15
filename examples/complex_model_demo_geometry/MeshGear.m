function meshobj=MeshGear(M,geom,mesh,specimen,meshobj,meshsizemin,meshsize)

  % geom=M.geom; mesh=M.mesh; meshobj=ModelWrapper(M,'specimen_mesh'); 
  %   or
  % geom=M.geom; mesh=M.mesh; meshobj=FindWrappedObject(M,'specimen_mesh');
  % specimen=FindWrappedObject(M,'specimen');
  
  % imported drawing is converted to a solid and swept 180 degrees in both CCW
  % and CW directions


  % 7/10/20 - Previously the disk was meshed with a swept mesh as much as
  % possible, with a free tetrahedral mesh for the region that included the
  % crack. For simplicity and in the interest of saving time this was changed
  % to a completely free tetrahedral mesh. Jared Taylor

  specimendomains=GetAutomaticSelectionEntities(M,geom,specimen,'dom');
  addprop(meshobj,'volume');
  meshobj.volume=ModelWrapper(M,[meshobj.tag '_volume']);
  BuildMeshFreeTet(M,geom,mesh,specimen,meshobj,meshsizemin,meshsize,@(M,geom,object) specimendomains);

  % Should probably separate out the size parameters for the different meshing steps

  % % This results in brick-like elements.
  % % CreateWrappedProperty(M,meshobj,'sweep',[meshobj.tag '_sweep'],mesh.node,'FreeTri');
  % CreateWrappedProperty(M,meshobj,'sweep',[meshobj.tag '_sweep'],mesh.node,'Sweep''FreeTri');
  % % meshobj.sweep.node.selection.geom(geom.tag,3); % select in 3-space (domain)

  % specimendomains=GetAutomaticSelectionEntities(M,geom,specimen,'dom');

  % % Determine which domain(s) make up the flange flat and which make up the
  % % remainder of the disk. The total volume of the flat is 48e-6 m^3 whereas
  % % the total volume of the remainder is 460e-6 m^3. We assume that all domains
  % % less than 100e-6 m^3 are part of the flange flat
 
  % [flangeflatdomains,bulkdomains]=FilterMatching(specimendomains, @(domainentity) DomainMeasureVolume(M,geom,domainentity) < 500e-6);

  % meshobj.sweep.node.selection.set(bulkdomains); % Domain selection?
  % % meshobj.sweep.node.set('sweeppath','circular');

  % % Comment for transducer, uncomment for welder
  % meshobj.sweep.node.set('facemethod','tri');

  % % Source face is the face created by the drawing import and conversion to
  % % solid.
  % % We get the automatic selection entities from specimen.crosssect.
  % % Half of these are the from the source face (say x > 0) and half are
  % % from the dest face (say x < 0). 
  % % So we get adjacent points
  % crosssect_boundaries=GetAutomaticSelectionEntities(M,geom,specimen.crosssect,'bnd');

  % [crosssect_source_boundaries,crosssect_dest_boundaries]=FilterMatching(crosssect_boundaries, @(crosssect_boundary) BoundaryEvaluateAtFirstVertexCoord(M,geom,crosssect_boundary,@(vertexcoord) vertexcoord(1) > 0.0));

  % % numbers of boundaries should match
  % assert(length(crosssect_source_boundaries)==length(crosssect_dest_boundaries))

  % % Filter the boundaries corresponding to the bulk of the cross section 
  % % from those corresponding to the flange flat. 
  % % 
  % % Do this by measuring the boundary's edge length. The Flange flat has
  % % an edge length of .0809 m whereas the bulk has an edge length of .6896
  % % 
  % % We switch at an edge length of 0.2
  % % 
  % % (We could also do this by intersecting with the boundary list of the
  % % previously determined flange domain)
  % [crosssect_source_bulk_boundaries,crosssect_source_flangeflat_boundaries]=FilterMatching(crosssect_source_boundaries, @(crosssect_boundary) BoundaryMeasureEdgeLength(M,geom,crosssect_boundary) > 0.2);
  % [crosssect_dest_bulk_boundaries,crosssect_dest_flangeflat_boundaries]=FilterMatching(crosssect_dest_boundaries, @(crosssect_boundary) BoundaryMeasureEdgeLength(M,geom,crosssect_boundary) > 0.2);
  % 
  % meshobj.sweep.node.selection('sourceface').set(crosssect_source_bulk_boundaries);
  % % End face is basically the same cross section, just 180 degrees rotated
  % % about the center.
  % meshobj.sweep.node.selection('targetface').set(crosssect_dest_bulk_boundaries);


  % % Create size node controlling size of meshing in sweep
  % CreateMeshSizeProperty(M,geom,mesh,meshobj.sweep,'size',[meshobj.sweep.tag '_size'],meshsizemin,meshsize,3,bulkdomains); % apply to our swept domain
  % 

  % % At the interface  between the swept
  % % region and the flange (squareish region with the bolt holes) , we need to convert to three sided facets
  % CreateWrappedProperty(M,meshobj,'convert',[meshobj.tag '_convert'],mesh.node,'Convert');
  % meshobj.convert.node.selection.geom(geom.tag,2); % 2d selection  
  % meshobj.convert.node.selection.set(FindContactBoundaries(M,geom,flangeflatdomains,bulkdomains));

  % % Triangular meshing of flange inner surface
  % CreateWrappedProperty(M,meshobj,'fl_inner',[meshobj.tag '_fl_inner'],mesh.node,'FreeTri');
  % % Identify inner surface by radius of all boundary vertices in x-y plane
  % % is less than xxx
  % flangeflatboundaries=mphgetadj(M.node,geom.tag,'boundary','domain',flangeflatdomains);
  % 
  % % The inner radius on the flange flat is at .1046 m, .1048 m and .1067 m
  % % The magnitude of the radii in the xy plain from the origin to all 
  % % of the vertices
  % % of a boundary B, collected into a matlab vector is less than .1080m 
  % %  sqrt(sum(bsxfun(@times,BoundaryEvaluateVertexCoords(M,geom,B),[1;1;0]).^2))
  % 
  % flangeflatinnerboundaries=FilterMatching(flangeflatboundaries,@(flangeflatboundary) all(sqrt(sum(bsxfun(@times,BoundaryEvaluateVertexCoords(M,geom,flangeflatboundary),[1;1;0]).^2)) < .1050));

  % 
  % meshobj.fl_inner.node.selection.set(flangeflatinnerboundaries);

  % % Create size node limiting size of meshing on inner boundary
  % CreateMeshSizeProperty(M,geom,mesh,meshobj.fl_inner,'size',[meshobj.fl_inner.tag '_size'],meshsizemin,meshsize,2,flangeflatinnerboundaries); % apply to our boundaries

  % % Triangular meshing of hole inner surface
  % CreateWrappedProperty(M,meshobj,'hole_inner',[meshobj.tag '_hole_inner'],mesh.node,'FreeTri');
  % % apply to boundaries of the extruded holes
  % meshobj.hole_inner.node.selection.set(GetAutomaticSelectionEntities(M,geom,specimen.holeextrude,'bnd'));

  % % Scaling for hole surface meshes.
  % meshobj.hole_inner.node.set('xscale','4');
  % meshobj.hole_inner.node.set('yscale','4');
  % meshobj.hole_inner.node.set('zscale','4');

  % % Volume meshing by FreeTet
  % % This handles remaining geometry of the specimen
  % % CreateWrappedProperty(M,meshobj,'volume',[meshobj.tag '_volume'],mesh.node,'FreeTet');
  % % meshobj.volume.node.selection.set(flangeflatdomains); 
  % addprop(meshobj,'volume');
  % meshobj.volume=ModelWrapper(M,[meshobj.tag '_volume']);
  % BuildMeshFreeTet(M,geom,mesh,specimen,meshobj,meshsizemin,meshsize,@(M,geom,object) flangeflatdomains);

