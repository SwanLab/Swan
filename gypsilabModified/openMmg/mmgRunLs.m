function mesh = mmgRunLs(mmg)
% Copyright (c) 2018-2019
% Matthieu Aussal, CMAP, Ecole Polytechnique
% Algiane Froehly, CARDAMOME, INRIA-SOFT
% LGPL Lesser General Public License v3.0.
% Remeshing using Mmg tools : https://www.mmgtools.org

% Current and mmg directory
here  = pwd;
there = which('mmg.m');
there = there(1:end-6);

% Move to mmg directory
cd(there)

% Check mesh (3D nodes, triangle or tetra, colours are int)
mesh = mmg.Mesh;
if (size(mesh.vtx,2) ~= 3) || (size(mesh.elt,2) < 3) || (size(mesh.elt,2) > 4)
    error('mmgRun.m : unavailable case, please use 3D vertices.')
end

% Check colours are int
if (norm(mesh.col-floor(mesh.col),'inf') > 1e-12) || (min(mesh.col)<0)
    error('mmgRun.m : only integer values for colours.')
end

% Delete .mesh and .sol files
deleteFiles()


%figure(1)
s.connec  = mesh.elt();
s.coord   = mesh.vtx(:,1:2);
mN1 = Mesh().create(s);
%mN1.plot();
msRelator = MasterSlaveRelator(s.coord);
masterSlave = msRelator.getRelation();
cell = CellNodesDescriptor(s.coord);
n1 = masterSlave(:);
n2 = cell.cornerNodes;
nT = [n1;n2];
mN1.computeEdges();
nodeInBoundaryEdges = mN1.edges.computeBoundaryEdges(nT);
edges = nodeInBoundaryEdges;
edges(:,3) = 0;
            

nodes = find(edges(:,3) ~= 10);

% Write original mesh in .mesh format
mmgMeshWrite('original.mesh',mesh,mmg.Req,nodes,edges);

% Write size map in .sol format
if ~isempty(mmg.Map)
    mmgSolWrite('original.sol',mmg.Map,is2d(mesh));
end

% Operating system
if ismac
    os = './mac_';
elseif ispc
    os = 'win_';
elseif isunix
    os = './uni_';
else
    disp('mmg.m : unavailable case')
end

% Element form
if (size(mesh.elt,2) == 3)
    if is2d(mesh)
        bin = 'mmg2d_O3';
    else
        bin = 'mmgs_O3';
    end
elseif (size(mesh.elt,2) == 4)
    bin = 'mmg3d_O3';
end

% Add extension (for windows only)
if ispc
    bin = [bin,'.exe '];
else
    bin = [bin,' '];
end

% Define temporary mesh files
file = '-in original.mesh -out refined.mesh ';

% Convert input to command
field = fieldnames(mmg);
opt   = '';
for i = 1:length(field)
    cmd = getfield(mmg,field{i});
    if ~isempty(cmd) && ~strcmp(field{i},'Mesh') && ~strcmp(field{i},'Req') && ~strcmp(field{i},'Map')
        opt = [opt,cmd,' '];
    end
end

if ~isempty(mmg.Map)
    cmd = ['-sol original.sol'];
    opt = [opt,cmd,' '];
    cmd = ['-ls 0'];
    opt = [opt,cmd,' '];
end


% Execute mmg binaries
command = [os,bin,file,opt,'-hausd 0.01 -hmin 0.001 -hmax 0.01 '];
%command = [os,bin,file,opt,'-hausd 0.005 -hmin 0.0005 -hmax 0.005 '];
system(command);

% Read refined mesh
[remesh,req,label,edges] = mmgMeshRead('refined.mesh');

% delete('refined.mesh');
% delete('original.sol');
% 
% it = remesh.col == 3;
% connec = remesh.elt(it,:);
% coord  = remesh.vtx(:,1:2);
% [s.coord,s.connec,edges] = computeUniqueCoordConnec(coord,connec,edges);
% 
% 
% 
% figure(1)
% clf
% mN2 = Mesh().create(s);
% mN2.plot()
% coord = mN2.coord;
% coord(:,3) = 0;
% remesh = msh(coord,mN2.connec);            
% 
% nodes = find(edges(:,3) ~= 10);
% mmgMeshWrite('original.mesh',remesh,[],nodes,edges);
% 
% command = [os,bin,file,' -hmin 0.001 -hmax 0.01 '];
% system(command);
% 
% [remesh,req,label,edges] = mmgMeshRead('refined.mesh');


% Replace required entites by original ones
if ~isempty(req)
    % Read original mesh for double/single troubles 
    origin = mmgMeshRead('original.mesh');
    
    % Vertex indices from required elements
    I = mesh.elt(mmg.Req,:);
    J = remesh.elt(req,:);

    % Unicity of vertices
    I = unique(I(:));
    J = unique(J(:));
    
    % Appairement using knnsearch with foat original mesh (to be improved)
    idx = knnsearch(origin.vtx(I,:),remesh.vtx(J,:),'K',1);
    
    % Replace remesh required vertices using original mesh vertices 
    remesh.vtx(J,:) = mesh.vtx(I(idx),:);
end

% Save output
mesh = remesh;

% Delete .mesh and .sol files
deleteFiles();

% Back to current directory
cd(here)
end

function [newCoord,newConnec,newEdges] = computeUniqueCoordConnec(coord,connec,edges)
allNodes = connec(:);
[uNodes,ind,ind2] = unique(allNodes,'rows','stable');

allCoords    = coord;
uniqueCoords = allCoords(uNodes,:);
newCoord     = uniqueCoords;

nnode = size(connec,2);
nCell = size(connec,1);
newConnec = reshape(ind2,nCell,nnode);

all(uNodes,1) = 1:length(uNodes);
newEdges(:,1) = all(edges(:,1));
newEdges(:,2) = all(edges(:,2));
newEdges(:,3) = edges(:,3);


end


function deleteFiles()
if exist('original.mesh','file')
    delete('original.mesh');
end
if exist('original.sol','file')
    delete('original.sol');
end
if exist('refined.mesh','file')
    delete('refined.mesh');
end
if exist('refined.sol','file')
    delete('refined.sol');
end
end
