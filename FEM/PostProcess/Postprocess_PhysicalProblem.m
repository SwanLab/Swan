classdef Postprocess_PhysicalProblem < Postprocess
    
    properties       
    end
    
    methods (Access = public)
        function obj = Postprocess_PhysicalProblem(results)
           obj.Values = results;
        end
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
        
        % Export Results 
        function ToGiDpost(obj,file_name,input,istep)
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
            obj.PrintVector(fid,ndim,'Displacements','U','Elastic Problem','Vector',1,'OnNodes','',obj.Values.d_u);
            obj.PrintTensor(fid,ndim,'Stress','S','Elastic Problem','Vector',1,'OnGaussPoints',gauss_points_name,obj.Values.stress);
            obj.PrintTensor(fid,ndim,'Strain','E','Elastic Problem','Vector',1,'OnGaussPoints',gauss_points_name,obj.Values.strain);
            fclose(fid);
        end       
    end
end

