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
                'testDualNestedInPrimalWithProjectedGradient';
                'testDualNestedInPrimalWithSlerp';                
                'testStressM1M2';
                'testM1M2';  
                'testVigdergauzMicroStructureWithStrain';                 
                'testVigdergauzMicroStructure';                                                
                'testMicro2';   
                'testMicro';                                                      
                'testBridge2';                
                'testGripping';   
                'testCantilever2';
                'testCantilever3';   
                'testCantilever';                
                'testBridge';                                         
                };
        end
        
    end
    
end
