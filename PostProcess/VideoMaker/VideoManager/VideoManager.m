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
            obj.createFolder();            
        end

        function makeVideo(obj,nIter)
            iterations = 0:nIter;
            obj.videoMaker.Set_up_make_video(iterations)
            obj.videoMaker.makeDesignVariableVideo();
        end
        
    end
    
    methods (Access = private)
        
        function createPaths(obj,cParams)
            fName = cParams.caseFileName;
            %obj.gidPath = 'C:\Program Files\GiD\GiD 13.0.4';% 
            obj.gidPath = '/opt/GiDx64/13.0.2/';
            obj.fileName = fName;
            obj.filePath = fullfile(pwd,'Output',fName);
        end
        
        function createVideoMaker(obj,s)
            cParams.type        = s.designVarType;
            cParams.dim         = s.pdim;
            cParams.fileName    = obj.fileName;
            cParams.filesFolder = obj.filePath;
            cParams.gidPath     = obj.gidPath;
            obj.videoMaker = VideoMaker.create(cParams);
        end
        
        function createFolder(obj)
            if ~exist(obj.filePath,'dir')
                mkdir(obj.filePath);
            end            
        end                        
        
    end
    
end

