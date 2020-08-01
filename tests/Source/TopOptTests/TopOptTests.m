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
                'testBridge2';                                              
                'testCantilever2';
                'testInteriorPerimeter';
                'testDualNestedInPrimalWithSlerp';
                'SimplAllTest3DExplicitVsImplicit';
                'SimplAllTest2DExplicitVsImplicit';
                'testMicro';
                'testDualNestedInPrimalWithProjectedGradient';
                'testCantilever';
                'testStressM1M2';
                'testM1M2';                             
                'testMicro2';
                'testVigdergauzMicroStructure';  
                'testSuperEllipseExponent';                                                
                'testVigdergauzMicroStructureWithStrain';
                'testAnalyticVsRegularizedPerimeter';
                'testGripping';  
                'testBridge';                                            
                };
        end
        
    end
    
end
