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
                'testTopOptDesignElemDensShapePrinting';                
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
