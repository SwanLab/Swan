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
                'testCantilever3';                                 
            'testMicro2';                  
                'testVigdergauzMicroStructureWithStrain';                 
                'testAnalyticVsRegularizedPerimeter';               
                'testVigdergauzMicroStructure';                                                                             
                'testGripping';   
                'testInteriorPerimeter';                
                'testDualNestedInPrimalWithSlerp';                  
                'testCantilever2';
                'SimplAllTest3DExplicitVsImplicit';                
                'SimplAllTest2DExplicitVsImplicit';
                'testMicro';                 
                'testDualNestedInPrimalWithProjectedGradient';
                'testCantilever';                
                'testStressM1M2';
                'testM1M2';               
                'testBridge';                                 
                 'testBridge2';                                                                                                
                };
        end
        
    end
    
end
