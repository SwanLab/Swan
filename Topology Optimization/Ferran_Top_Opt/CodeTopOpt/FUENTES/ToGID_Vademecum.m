function ToGID_Vademecum (file_name,istep,coordinatesa,conectivities,nnode)
% Construcci�n de malla postproceso

gtype = 'Triangle'; %gid type

nelem  = size(conectivities,1);           % Number of elements
npnod  = size(coordinatesa,1);                    % Number of nodes

msh_file = strcat(file_name,'_',num2str(istep),'.flavia.msh');

fid = fopen(msh_file,'w');
fprintf(fid,'### \n');
fprintf(fid,'# MAT_FEM  V.1.0 \n');
fprintf(fid,'# \n');

fprintf(fid,['MESH "WORKPIECE" dimension %3.0f   Elemtype %s   Nnode %2.0f \n \n'],2,gtype,nnode);
fprintf(fid,['coordinates \n']);


for i = 1 : npnod
    fprintf(fid,['%6.0f %12.5d %12.5d %12.5d \n'],i,coordinatesa(i,:));
end

fprintf(fid,['end coordinates \n \n']);
fprintf(fid,['elements \n']);

switch  gtype
    case 'Triangle'
        for i = 1 : nelem
            if (nnode==3)
            fprintf(fid,['%6.0f %6.0f %6.0f %6.0f  1 \n'],i,conectivities(i,:));
            end
        end
end

fprintf(fid,['end elements \n \n']);

status = fclose(fid);

end