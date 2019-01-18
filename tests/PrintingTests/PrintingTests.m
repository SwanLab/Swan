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
                  'testTopOptDesignElemDensShapePrinting';                                  
                  'testTopOptNonSelfAdjoint';                   
                  'testFEMPrinting';
                  'testTopOptDesignAndShapes';                                      
                  'testTopOptLevelSetGaussDensityPrinting';                                                                     
                  'testTopOptLevelSetPrinting';
                  'testTopOptDensityPrinting';   
                  'testTopOptGaussDensityPrinting';                                    
                };

        end
    end
    
end
