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
                'testMicro2';   
                'testMicro';                                      
                'testInteriorPerimeter';                                
                'testAnalyticVsRegularizedPerimeter';
                'testVigdergauzMicroStructure';                                
                'testVigdergauzMicroStructureWithStrain'; 
                'testDualNestedInPrimalWithProjectedGradient';
                'testDualNestedInPrimalWithSlerp';                
                'testStressM1M2';
                'testM1M2';  
                'testCantilever2';
                'testCantilever3';   
                'testCantilever';
                'testBridge2';                
                'testGripping';                                 
                'testBridge';                                                                
                };
        end
        
    end
    
end
