classdef Postprocess
    properties
        physical_problem
    end
    
    methods (Access = protected)
        function [nnode,ndim,pdim,gtype,etype,nelem,npnod,coordinates,conectivities] = getBasicParams(obj,input)
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
        function printTitle(obj,fid)
            fprintf(fid,'####################################################\n');
            fprintf(fid,'################# FEM-MAT-OO v.1.0 #################\n');
            fprintf(fid,'####################################################\n');
            fprintf(fid,'\n');
        end
        function PrintVector(obj,fid,ndim,nameres,indexName,problemType,result_type,istep,result_location,location_name,results)
            % Print Header ------------------------------------------------
            fprintf(fid,'\nResult "%s" "%s" %.0f %s %s "%s"\n',nameres,problemType,istep,result_type,result_location,location_name);
            switch ndim
                case 2
                    fprintf(fid,'ComponentNames  "%sx", "%sy"\n',indexName,indexName);
                case 3
                    fprintf(fid,'ComponentNames "%sx", "%sy", "%sz"\n',indexName,indexName,indexName);
                otherwise
                    error('Invalid value of parametre ndime.')
            end
            
            % Print Variables ---------------------------------------------
            fprintf(fid,'Values\n');
            for inode = 1:round(length(results)/ndim)
                fprintf(fid,'%6.0f ',inode);
                for idime = 1:ndim
                    fprintf(fid,'%12.5d ',results(ndim*(inode-1)+idime));
                end
                fprintf(fid,'\n');
            end
            fprintf(fid,'End Values\n');
        end
        function PrintTensor(obj,fid,ndim,nameres,indexName,problemType,result_type,istep,result_location,location_name,results)
            % Print Header ------------------------------------------------
            fprintf(fid,'\nResult "%s" "%s" %.0f %s %s "%s"\n',nameres,problemType,istep,result_type,result_location,location_name);
            switch ndim
                case 2
                    fprintf(fid,'ComponentNames  "%sx", "%sy", "%sxy", "%sz"\n',indexName,indexName,indexName,indexName);
                case 3
                    fprintf(fid,'ComponentNames "%sx", "%sy", "%sz", "%sxy", "%syz", "%sxz"\n',indexName,indexName,indexName,indexName);
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
        function PrintScalar(obj,fid,ndim,nameres,indexName,problemType,result_type,istep,result_location,location_name,results)
            % Print Header ------------------------------------------------
            fprintf(fid,'\nResult "%s" "%s" %.0f %s %s "%s"\n',nameres,problemType,istep,result_type,result_location,location_name);
            fprintf(fid,'ComponentNames  "%s"\n',indexName);
            
            % Print Variables ---------------------------------------------
            fprintf(fid,'Values\n');
            for inode = 1:length(results)
                fprintf(fid,'%6.0f ',inode);
                fprintf(fid,'%12.5d ',results(inode));
                fprintf(fid,'\n');
            end
            fprintf(fid,'End Values\n');
        end
        
        function  [fid,ndim,gauss_points_name] = write_header_res_file(obj,file_name,istep)
            
            res_file = fullfile('Output',strcat(file_name,'_',num2str(istep),'.flavia.res'));
            fid = fopen(res_file,'w');
            
            ngaus = obj.physical_problem.geometry.ngaus;
            posgp = obj.physical_problem.geometry.posgp';
            [~,ndim,~,~,etype,~,~,~,~] = obj.getBasicParams(obj.physical_problem);
            
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
            
        end
        
        
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
        
 
    end
    
    methods (Abstract, Access = protected)
        
        ToGiDpost(obj);
        
    end
    
    methods (Access = public)
        
        
        function  print(obj,file_name,istep)
            obj.ToGiD(file_name,obj.physical_problem,istep)
            obj.ToGiDpost(file_name,istep)
        end

        
        function file_list = create_file_list(obj,iterations_to_print,file_name,files_folder)
            file_list = [];
            for iter = 1:length(iterations_to_print)
                msh_file = fullfile(files_folder,strcat(file_name,'_',num2str(iterations_to_print(iter)),'.flavia.res'));
                file_list = [file_list, ' ',msh_file];
            end
            
        end
        
        function execute_tcl_files(obj,gidPath,file_tcl_name_with_path)
            system([fullfile(gidPath,'gid_offscreen'),' -t "source ',file_tcl_name_with_path,'"'])
            system(['rm ',file_tcl_name_with_path])
        end
        
        
    end
end

