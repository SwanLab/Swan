classdef GiDImageCapturer < handle

    properties (Access = private)
        resultsFile
        inputFileName
        outPutImageName
        gidPath
        swanPath
        pathTcl
        %outPutFolderPath
        outputImageName
    end

    methods (Access = public)

        function obj = GiDImageCapturer(cParams)
            obj.init(cParams);
            obj.createPathNames();
            obj.writeCallGiDTclFile();
        end

        function capture(obj)
            obj.captureImage();
            obj.cropImage();
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.resultsFile     = cParams.fileName;
            obj.outPutImageName = cParams.outPutImageName;
            obj.inputFileName   = cParams.inputFileName;
%             obj.swanPath = '/home/alex/git-repos/Swan/';
            obj.swanPath = '/home/ton/Github/Swan/';

%             obj.gidPath = '/home/alex/GiDx64/gid-15.0.3/';
            obj.gidPath = '/home/ton/GiDx64/gid-16.1.2d/';
        end

        function createPathNames(obj)
            obj.pathTcl = [obj.swanPath,'PostProcess/ImageCapturer/'];
            %obj.outPutFolderPath = [obj.swanPath,'Output/',obj.resultsFile,'/'];
            obj.outputImageName = [obj.outPutImageName];
        end

        function writeCallGiDTclFile(obj)
            tclFile = 'callGiDCapturer.tcl';
            obj.inputFileName = char(obj.inputFileName);
            %tclFileTocall = 'CaptureImage.tcl';
            tclFileTocall = 'CaptureImage3.tcl';
            % tclFileTocall = 'CaptureImageColor.tcl';
            %  tclFileTocall = 'CaptureSmoothImageColor.tcl';


            fid = fopen([obj.pathTcl,tclFile],'w+');
            fprintf(fid,['set path "',obj.pathTcl,'"\n']);
            fprintf(fid,['set tclFile "',tclFileTocall,'"\n']);
            fprintf(fid,['source $path$tclFile \n']);
            fprintf(fid,['set output ',obj.outputImageName,' \n']);
            fprintf(fid,['set inputFile ',obj.inputFileName,'\n']);
            fprintf(fid,['CaptureImage $inputFile $output \n']);
            fclose(fid);
        end

        function captureImage(obj)
            tclFile = [obj.pathTcl,'callGiDCapturer.tcl"'];
            command = [obj.gidPath,'gid_offscreen -offscreen -t "source ',tclFile];
            system(command);
        end

        function cropImage(obj)
            inputImage  = [' ',obj.outputImageName,'.png'];
            outPutImage = inputImage;
            % convert     = 'convert -crop 700x700+0+0 -gravity Center';
            convert     = 'convert -crop 442x442+0+0 -gravity Center';
            %convert     = 'convert -crop 1500x1500+0+0 -gravity Center';
            % convert     = 'convert -crop 600x300+0+0 -gravity Center';
            command = strcat(convert,' ',inputImage,' ',outPutImage);
            system(command);
        end

    end

end