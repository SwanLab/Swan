function [COOR,CONNECT,TypeElement1,CONNECTbound,TypeElementBOUND,MaterialType]=...
    ReadMeshFile(nameMSH,varargin)
%%% Leer .msh
% J.A. Hernandez
% % JAHO_B
% warning('JAHO_B')
% load('/home/joaquin/USO_COMUN_MATLAB/DATA/jaho.mat')
%dbstop('9')
if nargin == 0
    load('tmp.mat')
end
% INPUTS
READ_MATERIAL_COLUMN  = 0;
 %% EXTRACTING INPUTS
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
varginOR = varargin ;
FdnamesInputs = {'READ_MATERIAL_COLUMN'};
AuxDATA=WriteAuxFdNamesNEW(FdnamesInputs,varginOR);
for id = 1:length(AuxDATA);
    eval(AuxDATA{id});
end
%E_JAHO


%3) Open file
if strcmp(nameMSH(end-3:end),'.msh') == 0
    nameMSH = [nameMSH,'.msh']; 
end

fid=fopen(nameMSH,'r');

% First Read
%%%%%%%%%%%%
[OkFindRes RestLine] = ReadUntilToken(fid,'MESH'); % Find word "MESH"
% Data of the first mesh
[NameMesh1,ndime1,TypeElement1,nnode_elem1]=ObtInfMsh(RestLine);
%ndime1 = 3;
%nnode_elem1 = 4 ;
leido=fscanf(fid,'%s',1); %read "coordinates"
[OkFindRes ,num_nodos1] = ReadUntCount(fid,'end',ndime1+1) ;% counting nodes of mesh 1
[OkFindRes RestLine] = ReadUntilToken(fid,'Elements');
%leido=fscanf(fid,'%s',2); % read "coordinates" and "elements"
[OkFindRes ,num_elementos_1] = ReadUntCount(fid,'End',nnode_elem1+1+READ_MATERIAL_COLUMN); %% counting elements of mesh 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   boundary mesh
[OkFindRes RestLine] = ReadUntilToken(fid,'MESH'); % Find word "MESH"
if OkFindRes == 1 
% Data of the first mesh
[NameMesh1,ndime1,TypeElementBOUND,nnode_elemBOUND]=ObtInfMsh(RestLine);
[OkFindRes RestLine] = ReadUntilToken(fid,'Elements');

%leido=fscanf(fid,'%s',2); % read "coordinates" and "elements"
[OkFindRes ,num_elementosBOUND] = ReadUntCount(fid,'End',nnode_elemBOUND+1); %% counting elements of mesh 1
else
   % % there's no second mesh
   num_elementosBOUND = 0 ; 
   TypeElementBOUND = [] ; 
   nnode_elemBOUND = 0 ; 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Global Arrays
nmeshes = 1;
NameMeshes ={NameMesh1};
ndim       = ndime1;
nnode_elem = nnode_elem1;
num_elementos = num_elementos_1;
num_nodos = [num_nodos1];
%
%
% fclose(fid);
%
% % End of first lecture
%
% %%%%%%%% For storing information (structure)
%
% coordenadas = InitStrucField(NameMeshes,num_nodos,ndim+1);
% conexiones  = InitStrucField(NameMeshes,num_elementos,nnode_elem+1);
% material    = InitStrucField(NameMeshes,num_elementos,2*ones(nmeshes,1));

% Second lecture (storing)
fid=fopen(nameMSH,'r');
%% Lectura elementos y material(malla 1)
imesh = 1;
[OkFindRes RestLine] = ReadUntilToken(fid,'MESH');

% Coord
leido=fscanf(fid,'%s',1);  % read "coordinates"
COOR=leer_fichero_colum('',ndim(imesh)+1,0,num_nodos(imesh),0,0,fid);

COOR = COOR(:,2:end) ; 
if ~any(COOR(:,end)) 
 COOR = COOR(:,1:end-1);
end

% Conex
%dummy=fscanf(fid,'%s',3); % read "end coordinates" and "elements"
[OkFindRes RestLine] = ReadUntilToken(fid,'Elements');
CONNECT=leer_fichero_colum('',nnode_elem(imesh)+1+READ_MATERIAL_COLUMN,0,num_elementos(imesh),0,0,fid);
LIST_ELEMENTS = CONNECT(:,1) ; 
% Check if LIST_ELEMENTS is consecutive 
if any(diff(LIST_ELEMENTS)-1)
    warning('POSSIBLE MISTAKE IN READING THE MESH. SET "READ_MATERIAL_COLUMN" TO ZERO OR ASSIGN A MATERIAL' )
end
CONNECT = CONNECT(:,2:end) ; 

if READ_MATERIAL_COLUMN == 1
    MaterialType =  CONNECT(:,end) ; 
    CONNECT = CONNECT(:,1:end-1) ;
    
    
else
    MaterialType = [] ;
end

 

%
%% Boundary mesh  
imesh = 1;
if num_elementosBOUND > 0 
[OkFindRes RestLine] = ReadUntilToken(fid,'MESH');

% % Coord
% leido=fscanf(fid,'%s',1);  % read "coordinates"
% COOR=leer_fichero_colum('',ndim(imesh)+1,0,num_nodos(imesh),0,0,fid);

% Conex
%dummy=fscanf(fid,'%s',3); % read "end coordinates" and "elements"
[OkFindRes RestLine] = ReadUntilToken(fid,'Elements');
CONNECTbound=leer_fichero_colum('',nnode_elemBOUND+1,0,num_elementosBOUND,0,0,fid);
CONNECTbound = CONNECTbound(:,2:end) ; 
else
   CONNECTbound = [] ; 
   
end

% Jan-7-2019. Removing repeated elements 
if ~isempty(CONNECTbound)
CONNECTbound = CheckCONNECTbNONrepeated(CONNECTbound) ; 
end


fclose(fid);
