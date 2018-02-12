function ToGID (file_name,istep,coordinatesa,element,contactb,problembsc,nnode)
% Construcci�n de malla postproceso

etype = element.type;
ptype = problembsc.problemtype;
switch  etype
    case 'TRIANGLE'
        gtype = 'Triangle'; %gid type
    case 'QUAD'
        gtype = 'Quadrilateral';
    case 'HEXAHEDRA'
        gtype = 'Hexahedra';
end
nelem  = size(element.conectivities,1);           % Number of elements
npnod  = size(coordinatesa,1);                    % Number of nodes

msh_file = strcat(file_name,'_',num2str(istep),'.flavia.msh');

fid = fopen(msh_file,'w');
fprintf(fid,'### \n');
fprintf(fid,'# MAT_FEM  V.1.0 \n');
fprintf(fid,'# \n');

fprintf(fid,['MESH "WORKPIECE" dimension %3.0f   Elemtype %s   Nnode %2.0f \n \n'],2,gtype,nnode);
fprintf(fid,['coordinates \n']);
switch ptype
    case '2D'
        for i = 1 : npnod
            fprintf(fid,['%6.0f %12.5d %12.5d \n'],i,coordinatesa(i,:));
        end
    case '3D'
        for i = 1 : npnod
            fprintf(fid,['%6.0f %12.5d %12.5d %12.5d \n'],i,coordinatesa(i,:));
        end
end
fprintf(fid,['end coordinates \n \n']);
fprintf(fid,['elements \n']);

switch  gtype
    case 'Triangle'
%         for i = 1 : nelem
%             if (nnode==3)
%             fprintf(fid,['%6.0f %6.0f %6.0f %6.0f  1 \n'],i,element.conectivities(i,:));
%             end
%         end
        fprintf(fid,['%6.0f %6.0f %6.0f %6.0f  1 \n'],[1:nelem;element.conectivities']);
        
    case 'Quadrilateral'
%         for i = 1 : nelem
%             if (nnode==4)
%             fprintf(fid,['%6.0f %6.0f %6.0f %6.0f %6.0f  1 \n'],i,element.conectivities(i,:));
%             elseif (nnode==8)
%               fprintf(fid,['%6.0f %6.0f %6.0f %6.0f %6.0f  %6.0f %6.0f %6.0f %6.0f 1 \n'],i,element.conectivities(i,:));  
%             elseif (nnode==9)
%                fprintf(fid,['%6.0f %6.0f %6.0f %6.0f %6.0f %6.0f %6.0f %6.0f %6.0f %6.0f 1 \n'],i,element.conectivities(i,:)); 
%             end
%         end
        if (nnode==4)
            fprintf(fid,['%6.0f %6.0f %6.0f %6.0f %6.0f  1 \n'],[1:nelem;element.conectivities']);
        elseif (nnode==8)
            fprintf(fid,['%6.0f %6.0f %6.0f %6.0f %6.0f  %6.0f %6.0f %6.0f %6.0f 1 \n'],[1:nelem;element.conectivities']);
        elseif (nnode==9)
            fprintf(fid,['%6.0f %6.0f %6.0f %6.0f %6.0f %6.0f %6.0f %6.0f %6.0f %6.0f 1 \n'],[1:nelem;element.conectivities']);
        end

        
        
        
    case 'Hexahedra'
%         for i = 1 : nelem
%             fprintf(fid,['%6.0f %6.0f %6.0f %6.0f %6.0f %6.0f %6.0f %6.0f %6.0f  1 \n'],i,element.conectivities(i,:));
%         end
        
         fprintf(fid,['%6.0f %6.0f %6.0f %6.0f %6.0f %6.0f %6.0f %6.0f %6.0f  1 \n'],[1:nelem;element.conectivities']);
end

fprintf(fid,['end elements \n \n']);


status = fclose(fid);

end