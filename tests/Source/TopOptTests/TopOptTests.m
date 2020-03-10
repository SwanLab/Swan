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
                'testStressM1M2';
                'testM1M2';                  
                'testInteriorPerimeter';  
                'testDualNestedInPrimalWithProjectedGradient';
                'testDualNestedInPrimalWithSlerp';                
                'testVigdergauzMicroStructureWithStrain';                 
                'testVigdergauzMicroStructure';                                                              
                'testMicro';                                                      
                'testGripping';   
                'testCantilever2';
                'testCantilever';                
                'testAnalyticVsRegularizedPerimeter';               
                'testMicro2'; 
                'testBridge2';                
                'testCantilever3'; 
                'testBridge';                                                 
                };
        end
        
    end
    
end
