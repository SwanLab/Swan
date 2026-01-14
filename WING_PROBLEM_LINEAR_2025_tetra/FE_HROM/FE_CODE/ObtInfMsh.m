function [NameMesh1,ndime,TypeElement,nnode_elem]=ObtInfMsh(RestLine,varargin)
%        [NameMesh1,ndime,TypeElement,nnode_elem]=ObtInfMsh(RestLine)
% Read information about
%  MESH my_mesh2 dimension  2  Elemtype Triangle  Nnode  3
if nargin == 0 
    load('tmp.mat')
end

% DEFAULT INPUTS
COMFORM = 'NO';
%---------------
varginOR = varargin ;
FdnamesInputs = {'COMFORM'}; % MESH my_mesh dimension 2 Elemtype Triangle  Nnode 3 ... my_mesh field is missing
AuxDATA=WriteAuxFdNamesNEW(FdnamesInputs,varginOR);
for id = 1:length(AuxDATA);
    eval(AuxDATA{id});
end

switch COMFORM
    case 'NO'
NameMesh1 = [] ;
       % [NameMesh1 RestLine] = strtok(RestLine);
        [dummy,RestLine] = strtok(RestLine);
        [ndime,RestLine] = strtok(RestLine);
        if  ~isempty(ndime)
        ndime = str2num(ndime);
        else
           warning('ndime should be string ....Perhaps you forgot to assign material...') 
            ndime = str2num(ndime);
        end

        [dummy,RestLine] = strtok(RestLine);
        [TypeElement,RestLine] = strtok(RestLine);
        [dummy,RestLine] = strtok(RestLine);
        [nnode_elem,RestLine] = strtok(RestLine);
        %dbstop('17')
        if ~isempty(nnode_elem)
            nnode_elem = str2num(nnode_elem);
        else
            nnode_elem = [];
        end

    case 'YES'
        
        [NameMesh1 RestLine] = strtok(RestLine);
        [dummy,RestLine] = strtok(RestLine);
        [dummy,RestLine] = strtok(RestLine);
        [ndime,RestLine] = strtok(RestLine);
        ndime = str2num(ndime);

        [dummy,RestLine] = strtok(RestLine);
        [TypeElement,RestLine] = strtok(RestLine);
        [dummy,RestLine] = strtok(RestLine);
        [dummy,RestLine] = strtok(RestLine);
        [nnode_elem,RestLine] = strtok(RestLine);
        %dbstop('17')
        if ~isempty(nnode_elem)
            nnode_elem = str2num(nnode_elem);
        else
            nnode_elem = [];
        end


end