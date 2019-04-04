classdef VideoManager < handle
    
    properties (Access = private)
        videoMaker
        gidPath
        fileName
        filePath
    end
    
    methods (Access = public)
        
        function obj = VideoManager(settings,designVarType,pdim)
            obj.createVideoMaker(settings,designVarType,pdim);
            obj.gidPath = 'C:\Program Files\GiD\GiD 13.0.4';% 'C:\Program Files\GiD\GiD 13.0.3';
            obj.fileName = settings.case_file;
            obj.filePath = fullfile(pwd,'Output',settings.case_file);
        end
        
        function makeVideo(obj,nIter)
            iterations = 0:nIter;
            obj.videoMaker.Set_up_make_video(obj.gidPath,obj.fileName,obj.filePath,iterations)
            
            videoName = strcat('./Videos/Video_',obj.fileName,'_',int2str(nIter),'.gif');
            videoPath = fullfile(pwd,videoName);
            obj.videoMaker.Make_video_design_variable(videoPath)
        end
        
    end
    
    methods (Access = private)
        
        function createVideoMaker(obj,settings,designVarType,pdim)
            obj.videoMaker = VideoMakerTopOptFactory().create(settings.case_file,designVarType,pdim);
        end
        
    end
    
end

