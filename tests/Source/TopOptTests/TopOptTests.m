classdef TopOptTests < testRunner
    
    properties (Access = protected)
        FieldOfStudy = 'Topology Optimization'
        tests
    end
    
    methods (Access = public)
        
        function obj = TopOptTests()
            obj@testRunner();
        end
        
    end
    
    methods (Access = protected)
        
        function loadTests(obj)
            obj.tests = {... 
                'testCantilever2';
                 'testBridge2';                                                               
                'testCantilever3'; 
                'testDualNestedInPrimalWithProjectedGradient';
                'testDualNestedInPrimalWithSlerp';                  
                 'testMicro2';                 
                'testCantilever';                
                'testAnalyticVsRegularizedPerimeter';               
                'testInteriorPerimeter';                
                'testVigdergauzMicroStructureWithStrain';                 
                'testVigdergauzMicroStructure';                                                              
                'testMicro';                                                      
                'testGripping';   
                'testStressM1M2';
                'testM1M2';               
                'testBridge';                                 
                };
        end
        
    end
    
end
