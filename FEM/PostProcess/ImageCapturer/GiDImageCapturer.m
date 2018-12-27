classdef GiDImageCapturer
    
    properties (Access = private)
        gidPath
        resultsFile
        outputImageName
    end
    
    methods (Access = public)
        
        function obj = GiDImageCapturer(fileName,outPutName)            
            obj.init(fileName,outPutName)
        end
    end
    
    methods (Access = private)
        
        function init(obj,fileName,outPutName)
            obj.gidPath = '/opt/GiDx64/13.0.2/';    
            pathTcl = '/home/alex/git-repos/FEM-MAT-OO/FEM/PostProcess/ImageCapturer/';
            obj.resultsFile = fileName;
            obj.outputImageName = ['/home/alex/Dropbox/Amplificators/Images/',outPutName];            
            obj.writeCallGiDTclFile(pathTcl,fileName,obj.outputImageName);
            command = [obj.gidPath,'gid_offscreen -offscreen -t "source ',pathTcl,'callGiDCapturer.tcl"'];
            unix(command);
            obj.cropImage();
        end
        
        function writeCallGiDTclFile(obj,pathTcl,resultsFile,outputImageName)
                tclFile = 'callGiDCapturer.tcl';
                resultsFile = char(resultsFile);
                stlFileTocall = 'CaptureImage.tcl';
                fid = fopen([pathTcl,tclFile],'w+');
                fprintf(fid,['set path "',pathTcl,'"\n']);
                fprintf(fid,['set tclFile "',stlFileTocall,'"\n']);
                fprintf(fid,['source $path$tclFile \n']);
                fprintf(fid,['set output ',outputImageName,' \n']);
                fprintf(fid,['set inputFile ',resultsFile,'\n']);
                fprintf(fid,['CaptureImage $inputFile $output \n']);
                fclose(fid);
        end
         
        function cropImage(obj)
            name_file = [' ',obj.outputImageName,'.png'];
            command = strcat('convert -crop 500x500+0+0 -gravity Center ',name_file,' ',name_file);
            unix(command);
        end
        
            
        end
        
end
    

