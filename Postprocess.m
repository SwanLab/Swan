classdef Postprocess
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods(Static)
        % Mesh
        function ToGid(filename,iter,coord,connec,nnode,nelem,npnod)
            gtype = 'Linear';
            msh_file = strcat(filename,'_',num2str(iter),'.flavia.msh');
            fid = fopen(msh_file,'w');
            
            fprintf(fid,['MESH "WORKPIECE" dimension %2.0f  Elemtype %s   Nnode %2.0f \n \n'],3,gtype,nnode);
            fprintf(fid,['Coordinates \n']);
            for i = 1 : npnod
                fprintf(fid,['%6.0f %12.5d %12.5d %12.5d \n'],i,coord(i,:));
            end
            fprintf(fid,['End Coordinates \n \n']);
            fprintf(fid,['Elements \n']);
            
            fprintf(fid,['%6.0f %6.0f %6.0f \n'],[1:nelem;connec']);
            
            fprintf(fid,['End elements \n \n']);
            
            
            fclose(fid);
        end
        
        % Results
        function ToGidPost(filename,iter,ngaus,d_u)
            
            gtype = 'Linear';
            res_file = strcat(filename,'_',num2str(iter),'.flavia.res');
            fid = fopen(res_file,'w');
            job = 1;
            gid_write_headerpost(fid,gtype,ngaus,job)
            
            
            % DISPLACEMENT
            nameres = 'DISPLACEMENT';
            gid_write_vfield(fid,nameres,iter,d_u);
            
            fclose(fid);
        end
    end
    
end

