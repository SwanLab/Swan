classdef VideoMaker_Physical_Problem < VideoMaker
    
    
     properties 
         
         
     end
       
      methods (Access = public)
        function Make_video_stress(obj,output_video_name,component)
            post = Postprocess_PhysicalProblem;
            field2print = post.stress_name;
            componentfield = [post.stress_component,component];
            obj.Make_video_standard_field(field2print,componentfield,output_video_name)
        end
        
        function Make_video_strain(obj,output_video_name,component)
            post = Postprocess_PhysicalProblem;
            field2print = post.strain_name;
            componentfield = [post.strain_component,component];
            obj.Make_video_standard_field(field2print,componentfield,output_video_name)
        end
        
        function Make_video_displacement(obj,output_video_name,component)
            post = Postprocess_PhysicalProblem;
            field2print = post.displ_name;
            componentfield = [post.displ_component,component];
            obj.Make_video_standard_field(field2print,componentfield,output_video_name)
         end
        
      end
        
      methods (Access = private)
        function Make_video_standard_field(obj,field2print,componentfield,output_video_name_in)
            file_tcl_name = 'tcl_gid.tcl';
            file_list = obj.create_file_list(obj.iterations_to_print,obj.file_name,obj.files_folder);
            [output_video_name] = obj.replace_special_character(output_video_name_in);
            [file_list] = obj.replace_special_character(file_list);
            
            file_tcl_name_with_path = fullfile(obj.files_folder,file_tcl_name);
            file_path_in = fullfile(pwd,'FEM','PostProcess','Make_Video_stress.tcl');
            filepath = obj.replace_special_character(file_path_in);
                  
            fid = fopen(file_tcl_name_with_path,'w+');            

            fprintf(fid,'GiD_Process PostProcess \n');
            fprintf(fid,'%s\n',['set arg1 "',file_list,'"']);
            fprintf(fid,'%s\n',['set arg2 "',output_video_name,'"']);
            fprintf(fid,'%s\n',['set arg3 "',field2print,'"']);
            fprintf(fid,'%s\n',['set arg4 "',componentfield,'"']);  
                       
            fprintf(fid,'%s\n',['source "',filepath,'"']);
            fprintf(fid,'%s\n',['Make_Video_stress $arg1 $arg2 $arg3 $arg4']);
            fprintf(fid,'%s\n',['GiD_Process Mescape Quit']);
            fclose(fid);
            obj.executeTclFiles(obj.gidPath,file_tcl_name_with_path)
        end
        
      end

    
end