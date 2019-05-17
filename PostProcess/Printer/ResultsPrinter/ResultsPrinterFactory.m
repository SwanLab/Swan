classdef ResultsPrinterFactory < handle
    
    properties (Access = private)
        printer
    end
    
    methods (Access = public)
        function p = create(obj,resultCase,d)
            obj.createPrinter(resultCase,d);
            p = obj.printer;
        end
        
    end
    
    methods (Access = private)
        
        function createPrinter(obj,resultCase,d)
            switch resultCase
                case 'Elasticity'
                    p = ElasticityResultsPrinter(d);
                case 'ElasticityMicro'
                    p = ElasticityMicroResultsPrinter(d);
                case 'DensityGauss'
                    p = DensityGaussResultsPrinter(d);
                case 'LevelSet'
                    p = LevelSetResultsPrinter(d);
                case 'Density'
                    p = DensityResultsPrinter(d);
                case 'NumericalHomogenizer'
                    p = NumericalHomogenizerPrinter(d);
                case 'TopOptProblem'
                    p = TopOptResultsPrinter(d);
                case 'HomogenizedTensor'
                    p = HomogenizedTensorPrinter(d);
                case 'HomogenizedTensorStressBasis'
                    p = HomogenizedTensorStressBasisPrinter(d);
                case 'ElasticityMicroAndLevelSet'
                    p = ElasticityMicroDefinedByLevelSet(d);
                case {''}
                    p = [];
                case {'MicroParams'}
                    p = DensityResultsPrinter(d);
            end
            obj.printer = p;
        end
        
    end
    
end
