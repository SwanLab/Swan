classdef test2dQuad < testCheckStoredWithFemComputedVariable   
    properties (Access = protected)
        testName = 'test2d_quad';  
        variablesToStore = {'d_u'};
    end
       
    methods (Access = protected)
        
        function selectComputedVar(obj)
            obj.computedVar{1} = obj.fem.variables.d_u;            
        end        
        
    end
end

