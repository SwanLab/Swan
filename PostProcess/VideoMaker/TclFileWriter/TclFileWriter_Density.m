classdef TclFileWriter_Density < TclFileWriter
    
   methods (Access = public)
       
       function obj = TclFileWriter_Density(cParams)
           obj.tclTemplateName = 'Make_Video_density';                                           
           obj.init(cParams)           
       end
    
   end
   
     methods (Access = public)         

        function write(obj)
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