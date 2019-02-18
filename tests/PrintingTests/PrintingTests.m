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
                'testTopOptDesignElemDensShapePrinting';
                'testPrintingFreeFemFile';                
                'testPrintingInputFileForFem';
                'testNumericalHomogenizerPrinter';                                                                 
                'testFEMPrinting';
                'testTopOptLevelSetPrinting';
                'testTopOptDesignAndShapes';                
                'testTopOptNonSelfAdjoint';
                'testTopOptLevelSetGaussDensityPrinting';
                'testTopOptGaussDensityPrinting';
                'testTopOptDensityPrinting'; 
                };
            
        end
    end
    
end
