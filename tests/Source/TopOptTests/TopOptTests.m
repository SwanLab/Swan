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
               'testCantilever2';          
               'testCantilever';
               'testGripping';
               'testMicro';               
               'testMicro2';  
               'testDualNestedInPrimalWithProjectedGradient';                
               'testCantilever3';                
               'testDualNestedInPrimalWithSlerp';                           
               'testAnalyticVsRegularizedPerimeter';                
               'testSuperEllipseExponent';
               'testInteriorPerimeter';         
               'SimplAllTest3DExplicitVsImplicit';
               'SimplAllTest2DExplicitVsImplicit';                
               'testVigdergauzMicroStructureWithStrain';                                               
               'testVigdergauzMicroStructure';
                'testM1M2';                                                         
                'testStressM1M2';
                'testBridge2';               
                'testBridge';                              
                };
        end
        
    end
    
end
