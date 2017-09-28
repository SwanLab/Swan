classdef Postprocess
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods(Static)
        % Mesh
        function ToGid(file_name,input,istep)
            coordinates=input.mesh.coord;
            conectivities=input.mesh.connec;
            geometryType=input.mesh.geometryType;
            nnode=length(conectivities(1,:));
            ndim=input.mesh.ndim;
            etype = geometryType;
            ptype = '3D';
            switch  etype
                case 'Triangle'
                    gtype = 'Triangle'; %gid type
                case 'QUAD'
                    gtype = 'Quadrilateral';
                case 'Tetrahedra'
                    gtype = 'Tetrahedra';
                case 'HEXAHEDRA'
                    gtype = 'Hexahedra';
            end
            nelem  = size(conectivities,1);           % Number of elements
            npnod  = size(coordinates,1);                    % Number of nodes
            
            msh_file = strcat(file_name,'_',num2str(istep),'.flavia.msh');
            
            fid = fopen(msh_file,'w');
            fprintf(fid,'### \n');
            fprintf(fid,'# MAT_FEM  V.1.0 \n');
            fprintf(fid,'# \n');
            
            fprintf(fid,['MESH "WORKPIECE" dimension %3.0f   Elemtype %s   Nnode %2.0f \n \n'],ndim,gtype,nnode);
            fprintf(fid,['coordinates \n']);
            switch ptype
                case '2D'
                    for i = 1 : npnod
                        fprintf(fid,['%6.0f %12.5d %12.5d \n'],i,coordinates(i,:));
                    end
                case '3D'
                    for i = 1 : npnod
                        fprintf(fid,['%6.0f %12.5d %12.5d %12.5d \n'],i,coordinates(i,:));
                    end
            end
            fprintf(fid,['end coordinates \n \n']);
            fprintf(fid,['elements \n']);
            
            switch  gtype
                case 'Triangle'
                    %         for i = 1 : nelem
                    %             if (nnode==3)
                    %             fprintf(fid,['%6.0f %6.0f %6.0f %6.0f  1 \n'],i,conectivities(i,:));
                    %             end
                    %         end
                    fprintf(fid,['%6.0f %6.0f %6.0f %6.0f  1 \n'],[1:nelem;conectivities']);
                    
                case 'Quadrilateral'
                    %         for i = 1 : nelem
                    %             if (nnode==4)
                    %             fprintf(fid,['%6.0f %6.0f %6.0f %6.0f %6.0f  1 \n'],i,conectivities(i,:));
                    %             elseif (nnode==8)
                    %               fprintf(fid,['%6.0f %6.0f %6.0f %6.0f %6.0f  %6.0f %6.0f %6.0f %6.0f 1 \n'],i,conectivities(i,:));
                    %             elseif (nnode==9)
                    %                fprintf(fid,['%6.0f %6.0f %6.0f %6.0f %6.0f %6.0f %6.0f %6.0f %6.0f %6.0f 1 \n'],i,conectivities(i,:));
                    %             end
                    %         end
                    if (nnode==4)
                        fprintf(fid,['%6.0f %6.0f %6.0f %6.0f %6.0f  1 \n'],[1:nelem;conectivities']);
                    elseif (nnode==8)
                        fprintf(fid,['%6.0f %6.0f %6.0f %6.0f %6.0f  %6.0f %6.0f %6.0f %6.0f 1 \n'],[1:nelem;conectivities']);
                    elseif (nnode==9)
                        fprintf(fid,['%6.0f %6.0f %6.0f %6.0f %6.0f %6.0f %6.0f %6.0f %6.0f %6.0f 1 \n'],[1:nelem;conectivities']);
                    end
                case 'Tetrahedra'   
                    fprintf(fid,['%6.0f %6.0f %6.0f %6.0f %6.0f  1 \n'],[1:nelem;conectivities']);
                    
                    
                case 'Hexahedra'
                    %         for i = 1 : nelem
                    %             fprintf(fid,['%6.0f %6.0f %6.0f %6.0f %6.0f %6.0f %6.0f %6.0f %6.0f  1 \n'],i,conectivities(i,:));
                    %         end
                    
                    fprintf(fid,['%6.0f %6.0f %6.0f %6.0f %6.0f %6.0f %6.0f %6.0f %6.0f  1 \n'],[1:nelem;conectivities']);
            end
            
            fprintf(fid,['end elements \n \n']);
            
            
            status = fclose(fid);
            
        end
        
        % Results
        function ToGidPost(file_name,input,istep)
            geometryType=input.mesh.geometryType;
            ngaus = input.dim.ngaus;
            etype = geometryType;
            ndime=input.mesh.ndim; npnod=input.mesh.npnod; nnode=length(input.mesh.connec(1,:));
            nndof = input.dim.nunkn*npnod; 
            structural_values=input.variables;
            switch  etype
                case 'Triangle'
                    gtype = 'Triangle'; %gid type
                case 'QUAD'
                    gtype = 'Quadrilateral';
                case 'Tetrahedra'
                    gtype= 'Tetrahedra';
                case 'HEXAHEDRA'
                    gtype = 'Hexahedra';
            end
            
            
            % Escribe el fichero de resultados
            
            res_file =strcat(file_name,'_',num2str(istep),'.flavia.res');
            fid = fopen(res_file,'w');
            
            switch  etype
                case 'Triangle'
                    if nnode==3
                        idxgp = [1 2 3]; job=2;
                        gid_write_headerpost(fid,gtype,ngaus,job)
                    elseif nnode==6
                        idxgp = [];
                    end
                case 'QUAD'
                    if nnode==4
                        idxgp = [1 2 3 4 ]; job =3;
                        gid_write_headerpost(fid,gtype,ngaus,job)
                    elseif nnode==8
                        idxgp = [1 7 9 3 4 8 6 2 5]; job=1;
                        gid_write_headerpost(fid,gtype,ngaus,job)
                    end
                case 'Tetrahedra'
                    idxgp = [1 2 3 4 ]; job =3;
                    gid_write_headerpost(fid,gtype,ngaus,job)
                case 'HEXAHEDRA'
                    idxgp = [1 7 9 3 4 8 6 2 5]; job=1;
                    gid_write_headerpost(fid,gtype,ngaus,job)
            end
            
%             stres = structural_values.stress;
%             strain = structural_values.strain(:,1:nstre,:);
%             
%             nameres = ['Stress'];
%             postProcess.gid_write_gauss_tensorfield(fid,nameres,istep,stres,idxgp,ngaus,nstre,nelem);
%             
%             nameres = ['Strain'];
%             postProcess.gid_write_gauss_tensorfield(fid,nameres,istep,strain,idxgp,ngaus,nstre,nelem);
%            
%             
            nameres = 'DISPLACEMENT';
            tdisp(1:npnod,1)=structural_values.displacement(1:ndime:nndof);
            tdisp(1:npnod,2)=structural_values.displacement(2:ndime:nndof);
            if ndime==3
            tdisp(1:npnod,3)=structural_values.displacement(3:ndime:nndof);
            end
            
            gid_write_vfield(fid,nameres,istep,tdisp);
            
%             nameres = 'F_ext';
%             postProcess.gid_write_vfield(fid,nameres,istep,reshape(full(structural_values.fext),ndime,[])');
%             
           fclose(fid);

        end
    end
    
end

