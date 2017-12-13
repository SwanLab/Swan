classdef Postprocess
    properties
        Values
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
    end
    
    methods (Access = public)
        
        function Print_make_video_characteristic_function(obj,gidPath,file_name,files_folder,iterations_to_print,output_video_name)
        file_tcl_name = 'tcl_gid.tcl';
        field2print = 'LevelSet';
        componentfield = 'LS';
        file_list = obj.create_file_list(iterations_to_print,file_name,files_folder);
        min_value = -1e-32;
        max_value = 'Standard';
        file_tcl_name_with_path = fullfile(files_folder,file_tcl_name);
        fid = fopen(file_tcl_name_with_path,'w+');
        fprintf(fid,'GiD_Process PostProcess \n');
        fprintf(fid,['set arg1 "',file_list,'"\n']);
        fprintf(fid,['set arg2 "',output_video_name,'"\n']);
        fprintf(fid,['set arg3 "',field2print,'"\n']);
        fprintf(fid,['set arg4 "',componentfield,'"\n']);
        fprintf(fid,['set arg5 "',num2str(min_value),'"\n']);
        fprintf(fid,['set arg6 "',max_value,'"\n']);
        fprintf(fid,['source "',fullfile(pwd,'FEM','PostProcess','Make_Video_characteristic.tcl'),'"\n']);
        fprintf(fid,['Make_Video_characteristic $arg1 $arg2 $arg3 $arg4 $arg5 $arg6 \n']);
        fprintf(fid,['GiD_Process Mescape Quit']);
        fclose(fid);

        obj.execute_tcl_files(gidPath,file_tcl_name_with_path)
        end
        
        
        function Print_make_video_stress(obj,gidPath,file_name,files_folder,iterations_to_print,output_video_name)
            file_tcl_name = 'tcl_gid.tcl';
            field2print = 'Stress';
            componentfield = 'S';
            file_list = obj.create_file_list(iterations_to_print,file_name,files_folder);
            file_tcl_name_with_path = fullfile(files_folder,file_tcl_name);
            fid = fopen(file_tcl_name_with_path,'w+');
            fprintf(fid,'GiD_Process PostProcess \n');
            fprintf(fid,['set arg1 "',file_list,'"\n']);
            fprintf(fid,['set arg2 "',output_video_name,'"\n']);
            fprintf(fid,['set arg3 "',field2print,'"\n']);
            fprintf(fid,['set arg4 "',componentfield,'"\n']);
            fprintf(fid,['source "',fullfile(pwd,'FEM','PostProcess','Make_Video_stress.tcl'),'"\n']);
            fprintf(fid,['Make_Video_stress $arg1 $arg2 $arg3 $arg4 \n']);
            fprintf(fid,['GiD_Process Mescape Quit']);
            fclose(fid);

            obj.execute_tcl_files(gidPath,file_tcl_name_with_path)
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

