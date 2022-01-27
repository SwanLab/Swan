classdef DesignCapturerOfMaterialDesignExperiment < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        outputFigureName
        inputFile
        experimentPath
        pathInput
        iterations
    end
    
    properties (Access = private)
        
    end
    
    methods (Access = public)
        
        function obj = DesignCapturerOfMaterialDesignExperiment()
            obj.init()
            obj.captureAllImages();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj)
            path = '/home/alex/Dropbox/MaterialDesign/CC/Poisson';
            obj.outputFigureName = fullfile(path,'Poisson');
            obj.inputFile = 'HorizontalMaterialDesign';
            %obj.experimentPath = '/media/alex/MyPassport/MaterialDesign/CStar/LargeSquareWithPerimeter10_1';
            %obj.iterations = round(linspace(1,133,15));
            
         %   obj.experimentPath = '/media/alex/MyPassport/MaterialDesign/CStar/LargeCircleWithPerimeter10_12c';
         %   obj.iterations = round(linspace(1,67,15));
            
           % obj.experimentPath = '/media/alex/MyPassport/MaterialDesign/CStar/LargeCircleWithPerimeter10_1c';
           % obj.iterations = round(linspace(1,154,15));
         
            %obj.experimentPath = '/media/alex/MyPassport/MaterialDesign/CStar/HorizontalFromCircleWithPerimeter10_8_LargeCircle';
           % obj.iterations = round(linspace(1,27,9));
                     
            obj.experimentPath = '/media/alex/MyPassport/MaterialDesign/CStar/NegPoissonNoPerimeter5x6';
            obj.iterations = round(linspace(1,162,10));
            
            
            %obj.experimentPath = '/home/alex/git-repos/Swan/Output/HorizontalMaterialDesign';
            %obj.iterations = 0;
            
        end
        
        function captureAllImages(obj)
            nIter = length(obj.iterations);
            for iter = 1:nIter
                iterImage = obj.iterations(iter);
                obj.captureImage(iterImage);
            end
        end
                    
        function captureImage(obj,iter)
            i = iter;
            f = obj.inputFile;
            outPutNameWithIter = [obj.outputFigureName,'Iter',num2str(i)];
            inputFileName = fullfile(obj.experimentPath,[f,num2str(i),'.flavia.res']);
            s.fileName = f;
            s.outPutImageName = outPutNameWithIter;
            s.inputFileName   = inputFileName;
            imageCapturer     = GiDImageCapturer(s);
            imageCapturer.capture();
        end
        
    end
    
end