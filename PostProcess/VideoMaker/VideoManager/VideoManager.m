classdef VideoManager < handle
    
    properties (Access = private)
        videoMaker
        gidPath
        fileName
        filePath
    end
    
    methods (Access = public)
        
        function obj = VideoManager(cParams)
            obj.createPaths(cParams);
            obj.createVideoMaker(cParams);
        end
        
        function makeVideo(obj,nIter)
            iterations = 0:nIter;
            obj.videoMaker.Set_up_make_video(obj.gidPath,obj.fileName,obj.filePath,iterations)
            
            fileNameWithExtension = ['Video_',obj.fileName,'_',int2str(nIter),'.gif'];
            videoName = fullfile('Output',obj.fileName,fileNameWithExtension);
            videoPath = fullfile(pwd,videoName);
            
            if ~exist(obj.filePath,'dir')
                mkdir(obj.filePath);
            end
            
            obj.videoMaker.Make_video_design_variable(videoPath)
        end
        
    end
    
    methods (Access = private)
        
        function createVideoMaker(obj,cParams)
            type = cParams.designVarType;
            pdim = cParams.pdim;
            obj.videoMaker = VideoMakerTopOptFactory().create(obj.fileName,type,pdim);
        end
        
        function createPaths(obj,cParams)
            fName = cParams.caseFileName;
            %obj.gidPath = 'C:\Program Files\GiD\GiD 13.0.4';% 'C:\Program Files\GiD\GiD 13.0.3';
            obj.gidPath = '/opt/GiDx64/13.0.2/';% 'C:\Program Files\GiD\GiD 13.0.3';
            obj.fileName = fName;
            obj.filePath = fullfile(pwd,'Output',fName);
        end
        
    end
    
end

