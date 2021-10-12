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
                %% testTopOptCheckingDesignVariable
%                'testBridge2';
%                'testBridge';
%                'testCantilever2';
%                'testCantilever';
%                'testGripping';
%                'testMicro';
%                'testMicro2';
%                'testDualNestedInPrimalWithProjectedGradient';
%                'testCantilever3';
%                'testInteriorPerimeter';
%                'testDualNestedInPrimalWithSlerp';
                %%
%                'testAnalyticVsRegularizedPerimeter'; %tipusdiferent

                %%
%                'testSuperEllipseExponent';
                %%
%                'SimplAllTest3DExplicitVsImplicit';
%                'SimplAllTest2DExplicitVsImplicit';
                %%
%                'testVigdergauzMicroStructureWithStrain';
%                'testVigdergauzMicroStructure';
                %%
%                 'testM1M2'; % no funciona
%                 'testStressM1M2'; % no funciona
                };
        end
        
    end
    
end
