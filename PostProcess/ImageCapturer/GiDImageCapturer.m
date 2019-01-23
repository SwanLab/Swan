classdef GiDImageCapturer
    
    properties (Access = private)
        gidPath
        resultsFile
        outputImageName
    end
    
    methods (Access = public)
        
        function obj = GiDImageCapturer(fileName,outPutImageName,inputFileName)            
            obj.init(fileName,outPutImageName,inputFileName)
        end
    end
    
    methods (Access = private)
        
        function init(obj,fileName,outPutImageName,inputFileName)
            obj.gidPath = '/opt/GiDx64/13.0.2/';    
            pathTcl = '/home/alex/git-repos/FEM-MAT-OO/FEM/PostProcess/ImageCapturer/';
            obj.resultsFile = fileName;            
            obj.outputImageName = ['/home/alex/Dropbox/Amplificators/Images/',outPutImageName];            
            obj.writeCallGiDTclFile(pathTcl,inputFileName,obj.outputImageName);
            command = [obj.gidPath,'gid_offscreen -offscreen -t "source ',pathTcl,'callGiDCapturer.tcl"'];
            system(command);
            obj.cropImage();
        end
        
        function writeCallGiDTclFile(obj,pathTcl,inputFile,outputImageName)
                tclFile = 'callGiDCapturer.tcl';
                inputFile = char(inputFile);
                stlFileTocall = 'CaptureImage.tcl';
                fid = fopen([pathTcl,tclFile],'w+');
                fprintf(fid,['set path "',pathTcl,'"\n']);
                fprintf(fid,['set tclFile "',stlFileTocall,'"\n']);
                fprintf(fid,['source $path$tclFile \n']);
                fprintf(fid,['set output ',outputImageName,' \n']);
                fprintf(fid,['set inputFile ',inputFile,'\n']);
                fprintf(fid,['CaptureImage $inputFile $output \n']);
                fclose(fid);
        end
         
        function cropImage(obj)
            name_file = [' ',obj.outputImageName,'.png'];
            command = strcat('convert -crop 500x500+0+0 -gravity Center ',name_file,' ',name_file);
            system(command);
        end
        
            
        end
        
end
    

