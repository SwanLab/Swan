classdef PrintingTests < testRunner
    
    
    properties (Access = protected)
        FieldOfStudy = 'Printing'
        tests
    end
    
    methods (Access = public)
        function obj = PrintingTests()
            obj@testRunner();
        end
    end
    
    methods (Access = protected)
        function loadTests(obj)
            obj.tests = {...
                'testNumericalHomogenizerPrinter';                                                                 
                'testPrintingInputFileForFem';
                'testFEMPrinting';
                'testTopOptLevelSetPrinting';
                'testTopOptDesignAndShapes';                
                'testTopOptDesignElemDensShapePrinting';
                'testTopOptNonSelfAdjoint';
                'testTopOptLevelSetGaussDensityPrinting';
                'testTopOptGaussDensityPrinting';
                'testTopOptDensityPrinting'; 
                };
            
        end
    end
    
end
