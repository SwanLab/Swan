function gid_write_msh(fname,istep,coord,selem)
% 
% fname : file name
% istep : step number
% coord : mesh coordinates, coord(npnod,ndime)
% selem : structure containing the following info about the mesh
%         selem.matno = element material number
%         selem.conec = array of connectivities
%         selem.etype = element type, Triangle, Quadrilateral, ....

% check data
ok = check_data(selem);
if (ok==0) 
    return
end
% Initialization
nelem = size(selem.conec,1);  % Number of elements
nnode = size(selem.conec,2);  % Number of nodes per element
npnod = size(coord,1);        % Number of nodes
ndime = size(coord,2);        % Number of dimensions
    
conec = selem.conec;
matno = selem.matno;
eletyp = selem.etype;

% heading of the mesh file
msh_file = strcat(fname,'_',num2str(istep),'.flavia.msh');
fid = fopen(msh_file,'w');
fprintf(fid,'### \n');
fprintf(fid,'# ENG_AERO_COMP V.1.0 \n');
fprintf(fid,'# \n');
 
% nodal coordinates
fprintf(fid,['MESH "BODY" dimension %3.0f   Elemtype %s   Nnode %2.0f \n \n'],ndime,eletyp,nnode);
fprintf(fid,['coordinates \n']);
for i = 1 : npnod
    fprintf(fid,['%6.0f %12.5d %12.5d \n'],i,coord(i,:));
end
fprintf(fid,['end coordinates \n \n']);

% connectivities
fprintf(fid,['elements \n']);
for i = 1 : nelem
    fprintf(fid,['%6.0f %6.0f %6.0f %6.0f %6.0f \n'],i,conec(i,:),matno(i));
end
fprintf(fid,['end elements \n \n']);

% end of file
status = fclose(fid);
    
end