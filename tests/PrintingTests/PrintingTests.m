classdef PrintingTests < testRunner
    
    
    properties (Access = protected)
        FieldOfStudy = 'Printing Tests'
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
                'testTopOptNonSelfAdjoint';                                      
                'testTopOptDesignElemDensShapePrinting';                
                'testTopOptDesignAndShapes';         
                'testTopOptLevelSetGaussDensityPrinting';
                'testTopOptGaussDensityPrinting';
                'testTopOptDensityPrinting';
                'testFEMPrinting';
                'testTopOptLevelSetPrinting';               
                'testNumericalHomogenizerPrinter';                                                
                };
            
        end
    end
    
end
