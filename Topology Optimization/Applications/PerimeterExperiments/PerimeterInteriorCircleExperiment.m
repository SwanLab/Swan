classdef PerimeterInteriorCircleExperiment < handle
    
    methods (Access = public)
        
        function obj = PerimeterInteriorCircleExperiment()
            obj.computeGradientVariationWithRadius();            
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
            radius = 0.25;
            s.circleCase = 'interior';   
            s.levelSetParams.fracRadius = radius;
            s.levelSetParams.type = 'circle';
            L = 1;
            halfSide = L/2;
            s.levelSetParams.fracRadius = radius/halfSide;               
            s.curvature = 1/radius;
            s.circleCase = 'interior';
            s.nameCase = 'GradientCirclePerimeterExperiment';
            nameRoot   = 'SquareMacroTriangle';
            s.inputFiles = {nameRoot; [nameRoot,'Fine'];[nameRoot,'FineFine']};
            s.outputFolder = '/home/alex/git-repos/Perimeter/AllImages/';                 
        end         
        
    end
    
end