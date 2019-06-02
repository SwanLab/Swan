classdef  VideoMaker_TopOpt_density < VideoMaker_TopOpt    
    
    methods (Access = public)
        
        function obj = VideoMaker_TopOpt_density(cParams)
            obj.fieldName = 'Density';  
            obj.tclTemplateName = 'Make_Video_density';                        
            obj.init(cParams);
        end
        
    end
    
    methods (Access = protected)         

        function writeTclFileName(obj)
            fid = fopen(obj.tclFileName,'w+');            
            fprintf(fid,'GiD_Process PostProcess \n');
            fprintf(fid,'%s\n',['set arg1 "',obj.fileList,'"']);
            fprintf(fid,'%s\n',['set arg2 "',obj.videoFileName,'"']);
            fprintf(fid,'%s\n',['set arg3 "',obj.fieldName,'"']);
            fprintf(fid,'%s\n',['set arg4 "',obj.fieldName,'"']);
            fprintf(fid,'%s\n',['set arg5 "',obj.photoFileName,'"']);            
            fprintf(fid,'%s\n',['source "',obj.fullTclTemplateName,'"']);
            fprintf(fid,'%s\n',[obj.tclTemplateName,' $arg1 $arg2 $arg3 $arg4 $arg5']);
            fprintf(fid,'%s\n',['GiD_Process Mescape Quit']);
            fclose(fid);
        end               
        
    end
    
end
