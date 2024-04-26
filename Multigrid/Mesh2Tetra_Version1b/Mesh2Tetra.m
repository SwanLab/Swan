function T=Mesh2Tetra(V,F,options)
% This function MESH2TETRA converts a triangulated surface mesh into a 
% tetrahedron volume mesh. 
%
% Main advantage above existing constrained 3D Delaunay is that it 
% will never add new boundary points, (usefull for active appearance models) 
% Disadvantage, some highly non-convex surfaceshapes cannot be converted.
%
%   T=Mesh2Tetra(V,F,options);
%
% inputs,
%   V : Vertex List N x 3, with x,y,z positions
%   F : Face List M x 3, with vertex indices
%   options : Struct with options
%   options.verbose : if true, show information
%   options.checkinput : if true, check input mesh on errors
%
% outputs,
%   T : Tetrahedron List K x 4, with tetrahedron indices
%
%
% Note!, most functions are also available as c-code (much faster), 
%   run compile_c_files.m to compile the code
%  
% How the software works:
% - First, normal Delaunay is used to created a tetrahedron convexhull. 
%   Then  outside tetrahedrons and tetrahedrons intersecting the boundary 
%   mesh are removed. 
% - Second, New triangulated surface meshes are constructed for the space 
%   not yet filled by tetrahedrons. After which Delaunay
%   is done on the new boundary meshes.
% - Third, The remaining boundary which cannot be filled using Delaunay
%   constraints, is filled with a "Boundary collapse method". The Boundary
%   collapse method merges vertex neighbors, creating tetrahedrons while
%   making the surface mesh smaller (like a deflating balloon)
% - Fourth, It is possible that a part of the boundary mesh is left over which
%   cannot be filled with tetrahedrons. This is the case if there are no 4
%   vertices left who can see each other (like a non-convex polygon). In 
%   that case nearby Tetrahedrons are removed creating a new boundary 
%   mesh. And tetrahedron fitting with the boundary collapse methods is 
%   tried again (until succes, or a fixed amound of tries).
%
% Example,
%  compile_c_files.m
%
%  % Load the Training Data
%  load('example_jaw');
%
%  % Do normal delaunay tetrahedron generation
%  Tn = delaunayn(FV.vertices);
%
%  % Do constrained tetrahedron generation
%  Tc = Mesh2Tetra(FV.vertices,FV.faces);
%
%  figure,
%  subplot(2,2,1), hold on;
%   patch(FV,'facecolor',[1 1 0]);
%   axis equal; title('Surface Mesh');
%   plot3(FV.vertices(:,1),FV.vertices(:,2),FV.vertices(:,3),'r*');
%  subplot(2,2,2), hold on;
%   tetramesh(Tn,FV.vertices,'facecolor',[0 1 0]);
%   axis equal; title('Tetrahedrons After Delaunayn');
%  subplot(2,2,3), hold on;
%   tetramesh(Tc,FV.vertices,'facecolor',[0 0 1]);
%   axis equal; title('Tetrahedrons After Mesh2Tetra');
%
%
% Function is written by D.Kroon University of Twente (June 2010)

%% Process inputs
% Make all functions available
functionname='Mesh2Tetra.m';
functiondir=which(functionname);
functiondir=functiondir(1:end-length(functionname));
addpath([functiondir '/functions'])
addpath([functiondir '/functions/subfunctions'])
addpath([functiondir '/functions/mexfunctions'])


% Check input options
defaultoptions=struct('verbose',true,'checkinput',true);
if(~exist('options','var')), options=defaultoptions;
else
    tags = fieldnames(defaultoptions);
    for i=1:length(tags), if(~isfield(options,tags{i})), options.(tags{i})=defaultoptions.(tags{i}); end, end
    if(length(tags)~=length(fieldnames(options))),
        warning('Mesh2Tetra:unknownoption','unknown options found');
    end
end

%% Check Input Mesh 
if(options.checkinput)
	intersect=CheckInputMesh(V,F);

    % Solve Mesh intersections
    if(intersect(1))
        if(options.verbose)
            disp('Solve Mesh Intersections, by merging local faces');
        end
        [F,V]=solveInterSections(F,V);
    end
    if(intersect(4))
        disp('Swap Face orientation');
        F=[F(:,3) F(:,2) F(:,1)];
    end
end

VolumeMeshStart=CheckVolumeFaceMesh(V,F);
if(options.verbose)
    disp(['Volume inside boundary mesh is : ' num2str(VolumeMeshStart)]);
end

%% Do Delaunay Triangulation
if(options.verbose)
    disp('.');
    disp('Start Delaunay to Tetrahedrons');
end

[T F]=DelaunayInside3D(V,F,struct('verbose',options.verbose));

% Check size of Tetra Volume
VolumeTetra = CheckVolumeTetraMesh(V,T);
% Check size of Left Boundary Mesh
VolumeDelaunayFinish = CheckVolumeFaceMesh(V,F);

diff=(VolumeDelaunayFinish+VolumeTetra)-VolumeMeshStart;
if(abs(diff>1e-8))
    error('Mesh2Tetra:DelaunayInside3D',['Difference : ' num2str(diff)]);
end
        
if(options.verbose)
    disp('Finished Delaunay Triangulation');
    disp(['Volume boundary mesh left : ' num2str(VolumeDelaunayFinish)]);
    disp('.');
    disp('Start Boundary Collapse to Tetrahedrons');
end

%% Do Boundary Collapse Triangulation
T=BoundaryCollapse3D(V,F,T,struct('verbose',options.verbose,'checkinput',false));

% Check size of Tetra Volume
VolumeTetra = CheckVolumeTetraMesh(V,T);

diff=VolumeTetra-VolumeMeshStart;
if(abs(diff>1e-8))
    error('Mesh2Tetra:BoundaryCollapse3D',['Difference : ' num2str(diff)]);
end

