classdef Postprocess
    %Postprocess Summary of this class goes here
    %   Detailed explanation goes here
    
    %% !! CONSIDER WHETHER OR NOT THIS SHOULD BE O.O. !!
    
    properties
    end
    
    methods (Access = public)
        % Export Mesh
        function ToGiD(obj,file_name,input,istep)
            [nnode,ndim,pdim,gtype,etype,nelem,npnod,coordinates,conectivities] = obj.getBasicParams(input);
            
            msh_file = fullfile('Output',strcat(file_name,'_',num2str(istep),'.flavia.msh'));
            
            fid = fopen(msh_file,'w');
            obj.printTitle(fid);
            
            fprintf(fid,'MESH "WORKPIECE" dimension %3.0f   Elemtype %s   Nnode %2.0f \n \n',ndim,etype,nnode);
            fprintf(fid,'coordinates \n');
            switch pdim
                case '2D'
                    for i = 1:npnod
                        fprintf(fid,'%6.0f %12.5d %12.5d \n',i,coordinates(i,1:ndim));
                    end
                case '3D'
                    for i = 1:npnod
                        fprintf(fid,'%6.0f %12.5d %12.5d %12.5d \n',i,coordinates(i,:));
                    end
            end
            fprintf(fid,'end coordinates \n \n');
            
            fprintf(fid,'elements \n');
            switch  gtype
                case 'TRIANGLE'
                    switch nnode
                        case 3
                            fprintf(fid,'%6.0f %6.0f %6.0f %6.0f  1 \n',[1:nelem;conectivities']);
                        case 6
                            fprintf(fid,'%6.0f %6.0f %6.0f %6.0f %6.0f %6.0f %6.0f 1 \n',[1:nelem;conectivities']);
                        otherwise
                            error('Element type %s with %.0f nodes has not being implemented.',gytpe,nnode);
                    end
                case 'QUAD'
                    switch nnode
                        case 4
                            fprintf(fid,'%6.0f %6.0f %6.0f %6.0f %6.0f  1 \n',[1:nelem;conectivities']);
                        case 8
                            fprintf(fid,'%6.0f %6.0f %6.0f %6.0f %6.0f  %6.0f %6.0f %6.0f %6.0f 1 \n',[1:nelem;conectivities']);
                        case 9
                            fprintf(fid,'%6.0f %6.0f %6.0f %6.0f %6.0f %6.0f %6.0f %6.0f %6.0f %6.0f 1 \n',[1:nelem;conectivities']);
                        otherwise
                            error('Element type %s with %.0f nodes has not being implemented.',gytpe,nnode);
                    end
                case 'TETRAHEDRA'
                    switch nnode
                        case 4
                            fprintf(fid,'%6.0f %6.0f %6.0f %6.0f %6.0f  1 \n',[1:nelem;conectivities']);
                        otherwise
                            error('Element type %s with %.0f nodes has not being implemented.',gytpe,nnode);
                    end
                case 'HEXAHEDRA'
                    switch nnode
                        case 8
                            fprintf(fid,'%6.0f %6.0f %6.0f %6.0f %6.0f %6.0f %6.0f %6.0f %6.0f  1 \n',[1:nelem;conectivities']);
                        otherwise
                            error('Element type %s with %.0f nodes has not being implemented.',gytpe,nnode);
                    end
                otherwise
                    error('Element type %s has not being implemented.',gytpe,nnode);
            end
            
            fprintf(fid,'end elements\n\n');
            fclose(fid);
        end
        
        % Export Results from Physical Problem
        function ToGiDpost(obj,file_name,input,istep)
            % Write the file with the results
            res_file = fullfile('Output',strcat(file_name,'_',num2str(istep),'.flavia.res'));
            fid = fopen(res_file,'w');
            
            ngaus = input.geometry.ngaus;
            posgp = input.geometry.posgp';
            [~,ndim,~,~,etype,~,~,~,~] = obj.getBasicParams(input);
            results = input.variables;
            
            switch  etype
                case 'TRIANGLE'
                    etype = 'Triangle'; %gid type
                case 'QUAD'
                    etype = 'Quadrilateral';
                case 'TETRAHEDRA'
                    etype= 'Tetrahedra';
                case 'HEXAHEDRA'
                    etype = 'Hexahedra';
            end
            
            %% File Header
            fprintf(fid,'GiD Post Results File 1.0\n\n');
            obj.printTitle(fid);
            
            %% Print Gauss Points Header
            gauss_points_name = 'Guass up?';
            fprintf(fid,'GaussPoints "%s" Elemtype %s\n',gauss_points_name,etype);
            fprintf(fid,'Number of Gauss Points: %.0f\n',ngaus);
            fprintf(fid,'Nodes not included\n');
            fprintf(fid,'Natural Coordinates: given\n');
            for igaus = 1:ngaus
                for idime = 1:ndim
                    fprintf(fid,'%12.5d ',posgp(igaus,idime));
                end
                fprintf(fid,'\n');
            end
            fprintf(fid,'End GaussPoints\n');
            
            %% Print Results
            obj.PrintVector(fid,ndim,'Displacements','Vector',1,'OnNodes','',results.d_u);
            obj.PrintTensor(fid,ndim,'Stress','Vector',1,'OnGaussPoints',gauss_points_name,results.stress);
            obj.PrintTensor(fid,ndim,'Strain','Vector',1,'OnGaussPoints',gauss_points_name,results.strain);
            fclose(fid);
        end
        
        %Export Results from any input
        function ToGiDpostX(obj,file_name,x,input,istep)
            % Write the file with the results
            res_file = fullfile('Output',strcat(file_name,'_',num2str(istep),'.flavia.res'));
            fid = fopen(res_file,'w');
            
            ngaus = input.geometry.ngaus;
            posgp = input.geometry.posgp';
            [~,ndim,~,~,etype,~,~,~,~] = obj.getBasicParams(input);
            
            switch  etype
                case 'TRIANGLE'
                    etype = 'Triangle'; %gid type
                case 'QUAD'
                    etype = 'Quadrilateral';
                case 'TETRAHEDRA'
                    etype= 'Tetrahedra';
                case 'HEXAHEDRA'
                    etype = 'Hexahedra';
            end
            
            %% File Header
            fprintf(fid,'GiD Post Results File 1.0\n\n');
            obj.printTitle(fid);
            
            %% Print Gauss Points Header
            gauss_points_name = 'Guass up?';
            fprintf(fid,'GaussPoints "%s" Elemtype %s\n',gauss_points_name,etype);
            fprintf(fid,'Number of Gauss Points: %.0f\n',ngaus);
            fprintf(fid,'Nodes not included\n');
            fprintf(fid,'Natural Coordinates: given\n');
            for igaus = 1:ngaus
                for idime = 1:ndim
                    fprintf(fid,'%12.5d ',posgp(igaus,idime));
                end
                fprintf(fid,'\n');
            end
            fprintf(fid,'End GaussPoints\n');
            
            %% Print Results
            obj.PrintVector(fid,ndim,file_name,'Vector',1,'OnNodes','',x);
            fclose(fid);
        end
    end
    
    methods (Access = private, Static)
        function [nnode,ndim,pdim,gtype,etype,nelem,npnod,coordinates,conectivities] = getBasicParams(input)
            coordinates = input.mesh.coord;
            conectivities = input.mesh.connec;
            gtype = input.mesh.geometryType;
            nnode = length(conectivities(1,:));
            ndim = input.dim.ndim;
            pdim = input.mesh.pdim;
            
            switch  gtype %gid type
                case 'TRIANGLE'
                    etype = 'Triangle';
                case 'QUAD'
                    etype = 'Quadrilateral';
                case 'TETRAHEDRA'
                    etype = 'Tetrahedra';
                case 'HEXAHEDRA'
                    etype = 'Hexahedra';
            end
            nelem = size(conectivities,1); % Number of elements
            npnod = size(coordinates,1);   % Number of nodes
        end
        
        function printTitle(fid)
            fprintf(fid,'####################################################\n');
            fprintf(fid,'################# FEM-MAT-OO v.1.0 #################\n');
            fprintf(fid,'####################################################\n');
            fprintf(fid,'\n');
        end
        
        function PrintVector(fid,ndim,nameres,result_type,istep,result_location,location_name,results)
            % Print Header ------------------------------------------------
            fprintf(fid,'\nResult "%s" "Elastic Problem" %.0f %s %s "%s"\n',nameres,istep,result_type,result_location,location_name);
            switch nameres
                case 'Displacements'
                    index = 'U';
                case 'Level-set'
                    index = 'LS';
            end
            switch ndim
                case 2
                    fprintf(fid,'ComponentNames  "%sx", "%sy"\n',index,index);
                case 3
                    fprintf(fid,'ComponentNames "%sx", "%sy", "%sz"\n',index,index,index);
                otherwise
                    error('Invalid value of parametre ndime.')
            end
            
            % Print Variables ---------------------------------------------
            fprintf(fid,'Values\n');
            for inode = 1:round(length(results)/ndim)-1
                fprintf(fid,'%6.0f ',inode);
                for idime = 1:ndim
                    fprintf(fid,'%12.5d ',results(ndim*(inode-1)+idime));
                end
                fprintf(fid,'\n');
            end
            fprintf(fid,'End Values\n');
        end
        function PrintTensor(fid,ndim,nameres,result_type,istep,result_location,location_name,results)
            % Print Header ------------------------------------------------
            fprintf(fid,'\nResult "%s" "Elastic Problem" %.0f %s %s "%s"\n',nameres,istep,result_type,result_location,location_name);
            switch nameres
                case 'Strain'
                    index = 'E';
                case 'Stress'
                    index = 'S';
            end
            switch ndim
                case 2
                    fprintf(fid,'ComponentNames  "%sx", "%sy", "%sxy", "%sz"\n',index,index,index,index);
                case 3
                    fprintf(fid,'ComponentNames "%sx", "%sy", "%sz", "%sxy", "%syz", "%sxz"\n',index,index,index,index);
                otherwise
                    error('Invalid value of parametre ndime.')
            end
            
            % Print Variables ---------------------------------------------
            fprintf(fid,'Values\n');
            for ielem = 1:size(results,3)
                fprintf(fid,'%6.0f ',ielem);
                for igaus = 1:size(results,1)
                    for istre = 1:size(results,2)
                        fprintf(fid,'%12.5d ',results(igaus,istre,ielem));
                    end
                    fprintf(fid,'\n');
                end
            end
            fprintf(fid,'End Values\n');
        end
    end
end

