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
                %%
               'testVigdergauzMicroStructureWithStrain';
               'testVigdergauzMicroStructure';
                %%
%                 'testM1M2'; % no funciona
%                 'testStressM1M2'; % no funciona
                };
        end
        
    end
    
end
