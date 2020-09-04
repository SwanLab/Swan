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
            'testDualNestedInPrimalWithSlerp';
            'testMicro2';  
                'testMicro';                  
                'testCantilever2'; 
                'testCantilever3';                                
                'testBridge2';                   
                'testGripping';                  
                'testInteriorPerimeter';
                'SimplAllTest3DExplicitVsImplicit';
                'SimplAllTest2DExplicitVsImplicit';
                'testDualNestedInPrimalWithProjectedGradient';
                'testCantilever';
                'testStressM1M2';
                'testM1M2';                             
                'testVigdergauzMicroStructure';  
                'testSuperEllipseExponent';                                                
                'testVigdergauzMicroStructureWithStrain';
                'testAnalyticVsRegularizedPerimeter';                 
                'testBridge';                                                            
                };
        end
        
    end
    
end
