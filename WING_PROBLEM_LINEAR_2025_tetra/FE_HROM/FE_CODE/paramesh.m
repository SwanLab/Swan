function paramesh(NameFile,COOR,CONNECT,NAMEPROJ,MaterialType,TypeElement)

 elem_type = TypeElement;
 NNode=size(CONNECT,2);
 ndime = size(COOR,2) ;
 nnod = size(COOR,1) ;
 npe =  size(CONNECT,2) ;
 nElem = size(CONNECT,1) ;
 
fid = fopen('GIDPOST/PROTOTYPEMESH.vtk','wt');

fprintf(fid,'# vtk DataFile Version 2.0\n');
fprintf(fid,['volume: ',NAMEPROJ,'\n']);
fprintf(fid,'ASCII\n');
fprintf(fid,'DATASET UNSTRUCTURED_GRID\n');
fprintf(fid,'POINTS %i double\n',nnod);

% Coordinates 
format_xnod = [];
for k=1:ndime,
   format_xnod = [format_xnod, '%15.5e'];
end
   format_xnod=[ format_xnod, '\n'];

fprintf(fid,format_xnod,COOR(:,1:ndime)');
fprintf(fid,'\n');
% Elements
fprintf(fid,'CELLS %i %i\n',nElem,nElem*(npe+1));

format_icone = [];
for k=1:npe,
   format_icone=[format_icone, '%7i'];
end
   format_icone=['%10i ' format_icone, '\n'];

CONNECT = CONNECT - ones(nElem,npe);   


%fprintf(fid,format_icone,[npe*ones(1,nElem)',CONNECT, MaterialType]');
fprintf(fid,format_icone,[npe*ones(1,nElem)',CONNECT]');

CELLTYPE = 12*ones(nElem,1);
format_icone = ['%10i' , '\n'];
fprintf(fid,'\n');
fprintf(fid,'CELL_TYPES %i \n',nElem);
fprintf(fid,format_icone,CELLTYPE);
fprintf(fid,'\n');


fprintf(fid,'CELL_DATA %i \n',nElem);
fprintf(fid,'SCALARS cell_scalars int 1\n');
fprintf(fid,'LOOKUP_TABLE Materials\n');

fprintf(fid,format_icone,MaterialType');

fclose(fid);
end