classdef test3dHexahedra < testCheckStoredWithFemComputedVariable
    
    
    properties (Access = protected)
        testName = 'test3d_hexahedra';  
        variablesToStore = {'d_u'};
    end
       
    methods (Access = protected)
        
        function selectComputedVar(obj)
            obj.computedVar{1} = obj.fem.variables.d_u;            
        end
        
    end
    

end

