classdef GradientVariationWithRadiusExperiment < handle
       
    properties (Access = private)
       iMesh
       gExperiment
    end
    
    properties (Access = private)
       nameCase
       inputFiles
       outputFolder
       levelSetParams         
    end
    
    methods (Access = public)
        
        function obj = GradientVariationWithRadiusExperiment(s)
            obj.init(s);
        end
        
        function compute(obj)
            for im = 1:numel(obj.inputFiles)
                obj.iMesh = im;
                obj.createGradientVariationExperiment();
                obj.computeGradientVariationWithRadius();
            end
        end
    end
        
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.nameCase       = cParams.nameCase;
            obj.inputFiles     = cParams.inputFiles; 
            obj.outputFolder   = cParams.outputFolder;
            obj.levelSetParams = cParams.levelSetParams;            
        end
        
        function createGradientVariationExperiment(obj)
            s.inputFile = obj.inputFiles{obj.iMesh};
            s.iMesh     = obj.iMesh;
            s.levelSetParams = obj.levelSetParams;            
            g = GradientVariationExperiment(s);
            obj.gExperiment = g;
        end
        
        function computeGradientVariationWithRadius(obj) 
            s.mesh                 = obj.gExperiment.backgroundMesh;
            s.regularizedPerimeter = obj.gExperiment.regularizedPerimeter;
            s.inputFile            = obj.inputFiles{obj.iMesh};
            s.nameCase             = obj.nameCase;
            s.outputFolder         = obj.outputFolder;
            s.domainLength         = obj.gExperiment.domainLength();            
            gComputer = GradientVariationWithRadiusComputer(s);
            gComputer.compute();
        end                 
        
    end
    
end