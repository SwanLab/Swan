classdef Plotting3DCasesHorizontal < handle
    
    properties (Access = private)
        finalIter
        path
    end
    
    methods (Access = public)
        
        function obj = Plotting3DCasesHorizontal()
            close all
            obj.init();
            obj.takeImages();  
            obj.takeVideo();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj)
            %obj.path = '/home/alex/Desktop/CantileverTetra/';
            obj.path = '/home/alex/git-repos/Swan/Output/CantileverFlat/';
            %path = '/home/alex/git-repos/Swan/Output/CantileverTetra/';
            obj.finalIter = 385;
        end
        
        function takeImages(obj)
          u = obj.computeUnfittedMesh(obj.finalIter);
          s.unfittedMesh = u;  
          s.outDir       = obj.path;
          i = ImageTakerFor3DShape(s);
          i.takeImages();
        end
        
        function u = computeUnfittedMesh(obj,iteration)             
            fullPath = [obj.path,'/UnfittedMesh'];           
            filePath = [fullPath,num2str(iteration),'.mat'];
            uC = load(filePath);
            u = uC.u;              
        end
        
        function takeVideo(obj)
           f1 = obj.createFramesFromOptimization();
           f2 = obj.createFramesFromTurningCamera();
           f = [f1,f2];
           obj.createVideoFromFrames(f)
        end
        
        function fr = createFramesFromOptimization(obj)
            f = figure(2);
            f.WindowState = 'maximized';  
            for iter = 1:obj.finalIter                
                u = obj.computeUnfittedMesh(iter);
                figure(f)
                clf(f)                
                u.plot();
                view(60,30)
                title(['Iter = ',num2str(iter)])
                drawnow
                fr(iter) = getframe(gcf) ;
            end            
        end
        
        function fr = createFramesFromTurningCamera(obj)
            nIter = 120;
            for iter = 1:nIter                
                [a,z] = view();
                a = a + (1/nIter)*(360);
                view(a,z)
                drawnow
                fr(iter) = getframe(gcf) ;
            end            
        end        
        
        function createVideoFromFrames(obj,frames)
            writerObj = VideoWriter([obj.path,'/iterations.avi'],'Uncompressed AVI');
            writerObj.FrameRate = 10;
            secsPerImage = 1;
            open(writerObj);
            for i=1:length(frames)
                frame = frames(i) ;
                for v=1:secsPerImage(1)
                    writeVideo(writerObj, frame);
                end
            end
            close(writerObj);            
        end
        
    end
    

    
end


