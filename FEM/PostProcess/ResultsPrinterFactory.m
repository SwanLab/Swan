classdef ResultsPrinterFactory < handle
    
    methods (Access = public)
        
        function obj = ResultsPrinterFactory()
        end
    end
    
    methods (Access = public, Static)
        function p = create(resultCase)
            
            switch resultCase
                case 'NodalDensity'
                    p = DensityResultsPrinter();
                case 'NodalLevelSet'
                    p = LevelSetResultsPrinter();
                case 'GaussDensity'
                    p = DensityGaussResultsPrinter();
                case 'Elasticity'
                    p = ElasticityResultsPrinter();
                case {'SLERP','PROJECTED SLERP', 'HAMILTON-JACOBI'}
                    p = LevelSetResultsPrinter;      
                case {'PROJECTED GRADIENT', 'MMA', 'IPOPT'}
                    p = DensityResultsPrinter;                    
            end            
            
        end
        
    end
end