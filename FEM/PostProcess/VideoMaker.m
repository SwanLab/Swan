classdef VideoMaker < handle
    
    
    properties
        gidPath
        files_folder
        iterations_to_print
        file_name
    end
    
    methods (Access = public)
        
        
        function Set_up_make_video(obj,gidPath,file_name,files_folder,iterations_to_print)
            obj.gidPath = gidPath;
            obj.file_name = file_name;
            obj.files_folder = files_folder;
            obj.iterations_to_print = iterations_to_print;
        end
        
    end
    
    methods (Static, Access = protected)
        
        function file_list = create_file_list(iterations_to_print,file_name,files_folder)
            file_list = [];
            
            for iter = 1:length(iterations_to_print)
                msh_file = fullfile(files_folder,strcat(file_name,'_',num2str(iterations_to_print(iter)),'.flavia.res'));
                file_list = [file_list, ' ',msh_file];
            end
           file_list= replace(file_list,'\','\\\\');
        end
        
        function execute_tcl_files(gidPath,file_tcl_name_with_path)
            %system([fullfile(gidPath,'gid_offscreen'),' -t "source ',file_tcl_name_with_path,'"'])
            file_tcl_name_tcl= replace(file_tcl_name_with_path,'\','\\');
            system(['"',fullfile(gidPath,'gid_offscreen'),'"', ' -t ' ,'"source ',file_tcl_name_tcl,'"'])
            system(['DEL ',file_tcl_name_with_path])
        end
        
        function [output_string] = replace_special_character(input_string)
            
            output_string = replace(input_string,'\','\\\\');
            
            
        end
        
        
    end
    
    
    
    
    
end