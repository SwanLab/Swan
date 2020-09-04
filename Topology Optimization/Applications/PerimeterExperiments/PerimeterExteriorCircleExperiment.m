classdef PerimeterExteriorCircleExperiment < handle
    
    methods (Access = public)
        
        function obj = PerimeterExteriorCircleExperiment()
            obj.computeGradientVariationWithBoundary();
            obj.computeGradientSurfPerimeterComputer();
        end
        
    end
    
    methods (Access = private)
        
        function computeGradientVariationWithBoundary(obj)
            s = obj.createParams();
            g = GradientVariationWithBoundaryExperiment(s);            
            g.compute();            
        end
        
        function computeGradientVariationWithRadius(obj)
            s = obj.createParams(); 
            g = GradientVariationWithRadiusExperiment(s);
            g.compute();
        end        
        
        function computeGradientSurfPerimeterComputer(obj)
            s = obj.createParams(); 
            g = GradientSurfPerimeterComputer(s);            
            g.compute();
        end
        
    end
        
    methods (Access = private, Static)
        
        function s = createParams()
            radius = 2;
            s.curvature = 1/radius;
            s.levelSetParams.type = 'full';            
            s.circleCase = 'exterior';
            s.nameCase = 'GradientCirclePerimeterExperiment';
            nameRoot = 'DoubleCircleMacro';
            s.inputFiles = {nameRoot; [nameRoot,'Fine'];[nameRoot,'FineFine']};
            s.outputFolder = '/home/alex/git-repos/Perimeter/AllImages/';                 
        end        
        
    end
    
end