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
                'testInteriorPerimeter';                                
                'testAnalyticVsRegularizedPerimeter';
                'testVigdergauzMicroStructure';                                
                'testVigdergauzMicroStructureWithStrain'; 
                'testDualNestedInPrimalWithProjectedGradient';
                'testDualNestedInPrimalWithSlerp';
                'testStressM1M2';
                'testM1M2';  
                'testCantilever2';
                'testMicro';                
                'testCantilever3';   
                'testCantilever';
                'testBridge2';                
                'testGripping';                 
                'testMicro2';                
                'testBridge';                                                
                };
        end
        
    end
    
end
