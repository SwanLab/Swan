classdef  VideoMaker_TopOpt_levelSet3DXmirror < VideoMaker_TopOpt_levelSet
  properties
      
      
      
  end
  
  methods (Access = public)

        function Make_video_design_variable(obj,output_video_name)
            post = Postprocess_TopOpt_levelSet();
            field2print = post.levelSet_name;
            componentfield = post.levelSet_name_component;
            obj.Make_video_characteristic_function(field2print,componentfield,output_video_name)
        end
        
        
  end
        
  methods (Access = private)      
         function Make_video_characteristic_function(obj,field2print,componentfield,output_video_name_in)
            file_tcl_name = 'tcl_gid.tcl';
            file_list = obj.create_file_list(obj.iterations_to_print,obj.file_name,obj.files_folder);
            
            [output_video_name] = obj.replace_special_character(output_video_name_in);
            [output_photo]=obj.replace_special_character(strcat(obj.files_folder,'\',obj.file_name,'.png'));
            [file_list] = obj.replace_special_character(file_list);            
            
            min_value = -1e-32;
            
            file_tcl_name_with_path = fullfile(obj.files_folder,file_tcl_name);
            file_path_in = fullfile(pwd,'FEM','PostProcess','Make_Video_characteristic3DXmirror.tcl');
            filepath = obj.replace_special_character(file_path_in);
            
            fid = fopen(file_tcl_name_with_path,'w+');            

            fprintf(fid,'GiD_Process PostProcess \n');
            fprintf(fid,'%s\n',['set arg1 "',file_list,'"']);
            fprintf(fid,'%s\n',['set arg2 "',output_video_name,'"']);
            fprintf(fid,'%s\n',['set arg3 "',field2print,'"']);
            fprintf(fid,'%s\n',['set arg4 "',componentfield,'"']);
            fprintf(fid,'%s\n',['set arg5 "',num2str(min_value),'"']);
            fprintf(fid,'%s\n',['set arg6 "',output_photo,'"']);
            
            
            fprintf(fid,'%s\n',['source "',filepath,'"']);
            fprintf(fid,'%s\n',['Make_Video_characteristic $arg1 $arg2 $arg3 $arg4 $arg5 $arg6']);
            fprintf(fid,'%s\n',['GiD_Process Mescape Quit']);
            fclose(fid);
            obj.execute_tcl_files(obj.gidPath,file_tcl_name_with_path)
        end
      
        
  end
end
