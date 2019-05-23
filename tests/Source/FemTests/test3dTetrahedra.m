classdef test3dTetrahedra < testCheckStoredWithFemComputedVariable
    
    
    properties (Access = protected)
        testName = 'test3d_tetrahedra';  
        variablesToStore = {'d_u'};
    end
       
    methods (Access = protected)
        
        function selectComputedVar(obj)
            obj.computedVar{1} = obj.fem.variables.d_u;            
        end
        
    end
    

end

