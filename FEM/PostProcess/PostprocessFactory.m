classdef PostprocessFactory < handle
    
    
    
    
    methods (Access = public)
        
        function obj = PostprocessFactory()
        end
    end
    
    methods (Access = public, Static)
        function p = create(Postprocess)
            
            switch Postprocess
                case 'NodalDensity'
                    p = Postprocess_TopOpt_density();
                case 'NodalLevelSet'
                    p = Postprocess_TopOpt_levelSet();
                case 'GaussDensity'
                    p = PostprocessDensityInGaussPoints();
                case 'Elasticity'
                    p = Postprocess_PhysicalProblem();
            end            
            
        end
        
    end
end