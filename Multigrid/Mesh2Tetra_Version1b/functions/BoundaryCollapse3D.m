function T=BoundaryCollapse3D(V,F,T,Options)
%   This function performs 2D triangulation or 3D tessellation. Advantage
%   above delaunay is that only returns the faces inside the object
%
%   T=boundary_collapse_3d(V,F,T,options);
%       V : Mx3 the vertices, with x,y,z coordinates
%       F : Nx3 faces with vertice indices
%		(optional)
%		T : Ox4 Existing Tetrahedrons 
%   	options : Struct with options
%   	options.verbose : if true, show information
%   	options.checkinput : if true, check input mesh on errors
%
%   Example 3D,
%    % Load the Training Data
%     load('../example_sphere');
%
%     figure,
%     subplot(2,2,1), hold on;
%     patch(FV,'facecolor',[1 1 0]);
%     axis equal; title('input surfaces mesh');
%
%    % Do normal delaunay triangulation
%     T = delaunayn(FV.vertices);
%
%    % Do boundary collapse tesselation
%     T2=BoundaryCollapse3D(FV.vertices,FV.faces);
%
%     plot3(FV.vertices(:,1),FV.vertices(:,2),FV.vertices(:,3),'r*');
%     subplot(2,2,2), hold on;
%     tetramesh(T,FV.vertices,'facecolor',[0 1 0]);
%     axis equal; title('tetrahedrons after delaunayn');
%     subplot(2,2,3), hold on;
%     tetramesh(T2,FV.vertices,'facecolor',[0 0 1]);
%     axis equal; title('tetrahedrons after boundary_collapse');
%
%
%  Function is written by D.Kroon University of Twente (April 2010)

% Check input Options
defaultOptions=struct('verbose',true,'checkinput',true);
if(~exist('Options','var')), Options=defaultOptions;
else
    tags = fieldnames(defaultOptions);
    for i=1:length(tags), if(~isfield(Options,tags{i})), Options.(tags{i})=defaultOptions.(tags{i}); end, end
    if(length(tags)~=length(fieldnames(Options))),
        warning('Mesh2Tetra:unknownoption','unknown Options found');
    end
end


if(nargin>2)
	nT=size(T,1); T(nT+1:size(F,1),:)=0;
else
    nT=0; T=zeros([size(F,1) 4]);
    
end

nF=size(F,1);
if(Options.checkinput)
	CheckInputMesh(V,F(1:nF,:),T(1:nT,:));
end

Volume_Original=CheckVolumeFaceMesh(V,F)+CheckVolumeTetraMesh(V,T(1:nT,:));

if(Options.verbose), disp('Start Mesh to Tetrahedrons'); end
retry=0;
per=0;
Options.mode=0;
while(nF>0)
    nTold=nT;
    [V,F,nF,T,nT]=collapse_edge(V,F,nF,T,nT,Volume_Original,Options);
    if(nTold==nT)
        Options.mode=1;
        if(Options.verbose), disp('Some region left which cannot be divided in tetrahedrons, remove some neighbour tetrahedrons and try again'); end
        retry=retry+1;
        [V,F,nF,T,nT]=retry_remove_tetrahedrons(V,F,nF,T,nT);
        if(mod(retry,5)==0)
            [V,F,nF,T,nT]=retry_remove_tetrahedrons(V,F,nF,T,nT);
        end
        if(mod(retry,10)==0)
            [V,F,nF,T,nT]=retry_remove_tetrahedrons(V,F,nF,T,nT);
        end
        if(retry>25)
            error('Tetrahedron fitting inside mesh failed (because some random vertice selection is involved retry-ing may help)'); 
        end
    end
	if(Options.verbose)
		pern=round((1-(nF/size(F,1)))*100);
		if(pern~=per), disp(['Done : ' num2str(pern) ' %']); drawnow;per=pern; end
	end
end
T=T(1:nT,:);
inter=VolumeCheck(V,F,nF,T,nT,Volume_Original);
if(inter)
    error('Tetrahedron Volume differs from Face Volume')
end
if(Options.verbose), disp('Tetrahedrons succesfull fitted'); end


      
                       