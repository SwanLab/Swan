classdef VideoMaker_Physical_Problem < VideoMaker
    
    
     properties 
         
         
     end
       
      methods (Access = public)
        function Make_video_stress(obj,output_video_name)
            post = Postprocess_PhysicalProblem;
            field2print = post.stress_name;
            componentfield = post.stress_component;
            obj.Make_video_standard_field(field2print,componentfield,output_video_name)
        end
        
        function Make_video_strain(obj,output_video_name)
            post = Postprocess_PhysicalProblem;
            field2print = post.strain_name;
            componentfield = post.strain_component;
            obj.Make_video_standard_field(field2print,componentfield,output_video_name)
        end
        
        function Make_video_displacement(obj,output_video_name)
            post = Postprocess_PhysicalProblem;
            field2print = post.displ_name;
            componentfield = post.displ_component;
            obj.Make_video_standard_field(field2print,componentfield,output_video_name)
         end
        
      end
        
        
      methods (Access = private)
        function Make_video_standard_field(obj,field2print,componentfield,output_video_name)
            file_tcl_name = 'tcl_gid.tcl';
            file_list = obj.create_file_list(obj.iterations_to_print,obj.file_name,obj.files_folder);
            file_tcl_name_with_path = fullfile(obj.files_folder,file_tcl_name);
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
            obj.execute_tcl_files(obj.gidPath,file_tcl_name_with_path)
        end
        
      end
    
  
    
    
end