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
                case {'ShapeFunction'}
                    p = ShapesPrinter(d);
                case {'PerimeterPrinter'}
                    p = PerimeterPrinter(d);
                case {'VectorGauss'}
                    d.fieldName = d.name;
                    p = VectorGaussVariablePrinter(d);
                case {'ScalarGauss'}
                    d.fieldName = d.name;
                    p = ScalarGaussVariablePrinter(d);
                case {'ScalarNodal'}
                    d.fieldName = d.name;
                    p = ScalarNodalVariablePrinter(d);
                case {'DensityGauss'}
                    p = DensityGaussResultsPrinter(d);
                case {'LevelSet'}
                    p = LevelSetResultsPrinter(d);
                otherwise
                    p = [];
            end
            obj.printer = p;
        end
        
    end
    
end
