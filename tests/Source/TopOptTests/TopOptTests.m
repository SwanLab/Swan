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
                'testM1M2';   
                
                'testCantilever';
                'testDualNestedInPrimalWithProjectedGradient';                
                'testCantilever3';                
                'testDualNestedInPrimalWithSlerp';
                'testCantilever2';                
                'testBridge2';
                'testGripping';
                'testStressM1M2';
                'testMicro';                             
                'testAnalyticVsRegularizedPerimeter';                
                'testSuperEllipseExponent';
                'testInteriorPerimeter';
                'SimplAllTest3DExplicitVsImplicit';
                'SimplAllTest2DExplicitVsImplicit';                                
                'testBridge';
                'testVigdergauzMicroStructure';
                'testVigdergauzMicroStructureWithStrain';                                
                };
        end
        
    end
    
end
