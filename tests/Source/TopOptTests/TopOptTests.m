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
                'testBridge2';                                                               
                'testBridge';                 
                'testCantilever2';
                'testCantilever';                
                'testAnalyticVsRegularizedPerimeter';               
                'testMicro2'; 
                'testCantilever3'; 
                'testInteriorPerimeter';  
                'testDualNestedInPrimalWithProjectedGradient';
                'testDualNestedInPrimalWithSlerp';                
                'testVigdergauzMicroStructureWithStrain';                 
                'testVigdergauzMicroStructure';                                                              
                'testMicro';                                                      
                'testGripping';   
                'testStressM1M2';
                'testM1M2';                 
                };
        end
        
    end
    
end
