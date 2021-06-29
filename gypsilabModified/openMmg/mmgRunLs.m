function mesh = mmgRunLs(mmg)
% Copyright (c) 2018-2019
% Matthieu Aussal, CMAP, Ecole Polytechnique
% Algiane Froehly, CARDAMOME, INRIA-SOFT
% LGPL Lesser General Public License v3.0.
% Remeshing using Mmg tools : https://www.mmgtools.org

% Current and mmg directory
%here  = pwd;
there = which('mmg.m');
there = there(1:end-6);

oldMeshFileName = mmg.oldMeshFileName;
newMeshFileName = mmg.newMeshFileName;
levelSetOldFileName = mmg.levelSetOldFileName;
levelSetNewFileName = mmg.levelSetNewFileName;

% Move to mmg directory
%cd(there)

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
deleteFiles(oldMeshFileName,newMeshFileName,levelSetOldFileName,levelSetNewFileName)


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
mmgMeshWrite(oldMeshFileName,mesh,mmg.Req,nodes,edges);

% Write size map in .sol format
if ~isempty(mmg.Map)
    mmgSolWrite(levelSetOldFileName,mmg.Map,is2d(mesh));
end

% Operating system
if ismac
    os = './mac_';
elseif ispc
    os = 'win_';
elseif isunix
    os = '/./uni_';
else
    disp('mmg.m : unavailable case')
end
os = [there,os];

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
file = ['-in ',oldMeshFileName,' -out ',newMeshFileName,' '];

% Convert input to command
field = fieldnames(mmg);
opt   = '';
for i = 1:length(field)
    cmd = getfield(mmg,field{i});
    if ~isempty(cmd) && ~strcmp(field{i},'Mesh') && ~strcmp(field{i},'Req')...
       && ~strcmp(field{i},'Map') ...
       && ~strcmp(field{i},'oldMeshFileName') && ~strcmp(field{i},'newMeshFileName') ...
       && ~strcmp(field{i},'levelSetOldFileName') && ~strcmp(field{i},'levelSetNewFileName')
       %&& ~strcmp(field{i},'original') && ~strcmp(field{i},'refined') ...
        opt = [opt,cmd,' '];
    end
end

if ~isempty(mmg.Map)
    cmd = ['-sol ',levelSetOldFileName];
    opt = [opt,cmd,' '];
    cmd = ['-ls 0'];
    opt = [opt,cmd,' '];
end


% Execute mmg binaries
command = [os,bin,file,opt];
[~,~] = system(command);

% Read refined mesh
[remesh,req,label,edges] = mmgMeshRead(newMeshFileName);

% Replace required entites by original ones
if ~isempty(req)
    % Read original mesh for double/single troubles 
    origin = mmgMeshRead(oldMeshFileName);
    
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
deleteFiles(oldMeshFileName,newMeshFileName,levelSetOldFileName,levelSetNewFileName);

% Back to current directory
%cd(here)
end


function deleteFiles(oldMeshFileName,newMeshFileName,levelSetOldFileName,levelSetNewFileName)
if exist(oldMeshFileName,'file')
    delete(oldMeshFileName);
end
if exist(newMeshFileName,'file')
    delete(newMeshFileName);
end
if exist(levelSetOldFileName,'file')
    delete(levelSetOldFileName);
end
if exist(levelSetNewFileName,'file')
    delete(levelSetNewFileName);
end

end
