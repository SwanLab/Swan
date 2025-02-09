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
               'testTopOptDesignAndShapes';                
               'testTopOptNonSelfAdjoint';
               'testTopOptLevelSetGaussDensityPrinting';
               'testTopOptGaussDensityPrinting';
                'testTopOptDensityPrinting'; 
               'testPrintingInputFileForFem';
                'testNumericalHomogenizerPrinter';                                                                                 
               'testTopOptDesignElemDensShapePrinting';
                'testPrintingFreeFemFile';                
                'testFEMPrinting';
                'testTopOptLevelSetPrinting';                
                };
            
        end
    end
    
end
