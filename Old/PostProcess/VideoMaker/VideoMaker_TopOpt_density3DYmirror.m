classdef  VideoMaker_TopOpt_density3DYmirror < VideoMaker_TopOpt
  properties
      
      
      
  end
  
  methods (Access = public)

        function Make_video_design_variable(obj,output_video_name)
            post = Postprocess_TopOpt_density();
            field2print = post.density_name;
            componentfield = post.density_name_component;
            obj.Make_video_density_field(field2print,componentfield,output_video_name)
        end
        
        
        function Make_video_design_variable_reg(obj,output_video_name)
            post = Postprocess_TopOpt_density();
            field2print = post.density_name_reg;
            componentfield = post.density_name_component_reg;
            obj.Make_video_density_field(field2print,componentfield,output_video_name)
        end
        
        
  end
        
  methods (Access = private)
        function Make_video_density_field(obj,field2print,componentfield,output_video_name)
            file_tcl_name = 'tcl_gid.tcl';
            file_list = obj.create_file_list(obj.iterations_to_print,obj.file_name,obj.files_folder);
            [output_video_name] = obj.replace_special_character(output_video_name);
            [output_photo]=obj.replace_special_character(strcat(obj.files_folder,'\',obj.file_name,'.png'));
            [file_list] = obj.replace_special_character(file_list);
            
            file_tcl_name_with_path = fullfile(obj.files_folder,file_tcl_name);
            file_path_in = fullfile(pwd,'FEM','PostProcess','Make_Video_density3DYmirror.tcl');
            filepath = obj.replace_special_character(file_path_in);
                  
            fid = fopen(file_tcl_name_with_path,'w+');            

            fprintf(fid,'GiD_Process PostProcess \n');
            fprintf(fid,'%s\n',['set arg1 "',file_list,'"']);
            fprintf(fid,'%s\n',['set arg2 "',output_video_name,'"']);
            fprintf(fid,'%s\n',['set arg3 "',field2print,'"']);
            fprintf(fid,'%s\n',['set arg4 "',componentfield,'"']); 
            fprintf(fid,'%s\n',['set arg5 "',output_photo,'"']);
                       
            fprintf(fid,'%s\n',['source "',filepath,'"']);
            fprintf(fid,'%s\n',['Make_Video_density $arg1 $arg2 $arg3 $arg4 $arg5']);
            fprintf(fid,'%s\n',['GiD_Process Mescape Quit']);
            fclose(fid);
            obj.execute_tcl_files(obj.gidPath,file_tcl_name_with_path)
        end
        
    end
end
