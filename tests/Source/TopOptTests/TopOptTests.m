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
                'testAnalyticVsRegularizedPerimeter';               
                'testMicro2'; 
                'testBridge2';                
                'testCantilever3';   
                'testInteriorPerimeter';                                                
                'testDualNestedInPrimalWithProjectedGradient';
                'testDualNestedInPrimalWithSlerp';                
                'testStressM1M2';
                'testM1M2';  
                'testVigdergauzMicroStructureWithStrain';                 
                'testVigdergauzMicroStructure';                                                              
                'testMicro';                                                      
                'testGripping';   
                'testCantilever2';
                'testCantilever';                
                'testBridge';                                         
                };
        end
        
    end
    
end
