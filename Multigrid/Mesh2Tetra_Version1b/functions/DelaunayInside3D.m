function [T F]=DelaunayInside3D(V,F,Options)
%
%  [T Fr]=DelaunayInside3D(V,F,Options)
%
%  inputs,
%    V : Vertex List N x 3, with x,y,z positions
%    F : Face List M x 3, with vertex indices
%    options : Struct with options
%    options.verbose : if true, show information
%
% outputs,
%    T : Tetrahedron List K x 4, with tetrahedron indices
%	 Fr : Remaining Faces List O x 3, with vertex indices
%
%  Function is written by D.Kroon University of Twente (April 2010)

% Check input Options
defaultOptions=struct('verbose',true,'itt',0);
if(~exist('Options','var')), Options=defaultOptions;
else
    tags = fieldnames(defaultOptions);
    for i=1:length(tags), if(~isfield(Options,tags{i})), Options.(tags{i})=defaultOptions.(tags{i}); end, end
    if(length(tags)~=length(fieldnames(Options))),
        warning('Mesh2Tetra:unknownoption','unknown Options found');
    end
end

os(1:1+Options.itt)='-';

% Triangulate the points inside the mesh
[nObjects Object]=ReturnSepparateFaceObjects(F);
if(Options.verbose), disp([os 'Number of objects : ' num2str(nObjects)]); end
T=[];
if(Options.verbose), disp([os 'Starting Volume : ' num2str(CheckVolumeFaceMesh(V,F)) ]); end
for i=1:nObjects,
    if(Options.verbose), disp([os 'Processing object : ' num2str(i)]); end
    F2 = Object{i};
    if(Options.verbose), disp([os 'Number of Faces : ' num2str(length(F2))]); end
    
    % Make a list which contains only the inside mesh points
    % and update the face list with the new vertex-id's
    [V2,F2,ID2] = InsidePoints3D(V,F2);
    
    %intersect=CheckMeshInterSections(V2,F2);
    %if(intersect),keyboard; end
        
    VolumeMeshStart = CheckVolumeFaceMesh(V2,F2);
    if(Options.verbose), disp([os 'Object Volume : ' num2str( VolumeMeshStart) ]); end
    
    % Do delaunay to get Faces
    T2 = delaunayn(V2);
    
    % Remove outside and other not valid Faces
    T2=RemoveInvalidTetrahedrons(T2,V2,F2);
       
    VolumeTetra = CheckVolumeTetraMesh(V2,T2);
    
    % Get the local not finished contours
    F2=GetRemainingFaces(T2,F2,V2);
    VolumeMeshFinish = CheckVolumeFaceMesh(V2,F2);
    
    if(Options.verbose), disp([os 'Remaining Face Volume : ' num2str(VolumeMeshFinish) ]); end
    diff=(VolumeMeshFinish+VolumeTetra)-VolumeMeshStart;
    if(abs(diff>1e-8))
        if(Options.verbose), disp([os 'Difference : ' num2str(diff)]); end
		T2=[];
    else
        intersect=CheckMeshInterSections(V2,F2);
        if(intersect), T2=[]; end
    end
    
    % Triangulate the remainig contours
    if((~isempty(F2))&&(~isempty(T2))), 
        if(Options.verbose), disp([os 'Remaining Number of Faces : ' num2str(length(F2))]); end
        Options2=Options; Options2.itt=Options2.itt+1;
        T3=DelaunayInside3D(V2,F2,Options2);
		T2=[T2;T3];
    end
        
    % Tetrahedrons ID's back to old list
    for j=1:size(T2,1)
        T2(j,1)=ID2(T2(j,1)); T2(j,2)=ID2(T2(j,2)); 
        T2(j,3)=ID2(T2(j,3)); T2(j,4)=ID2(T2(j,4));
    end

    % Add to total list
    T=[T;T2];
end
F=GetRemainingFaces(T,F,V);
intersect=CheckMeshInterSections(V2,F2);
if(intersect),
    warning('DelaunayInside3D:Process','Mesh Intersections Detected'); 
end
       
