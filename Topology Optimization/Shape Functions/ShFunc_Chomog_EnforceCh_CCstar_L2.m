classdef ShFunc_Chomog_EnforceCh_CCstar_L2 < ShFunc_Chomog_EnforceCh
  
    properties (Access = private)
        initial_value
        epsilon
    end
    
    methods (Access = public) 
        
        function obj=ShFunc_Chomog_EnforceCh_CCstar_L2(cParams)
            obj.initChomog(cParams);
            obj.computeChTarget(cParams.ChTarget);
            obj.epsilon = cParams.targetParameters.epsilon_isotropy;
        end
        
        function computeFunction(obj)
            %obj.computePhysicalData();
            obj.computeFunctionAndGradient();     
         %   obj.normalizeFunction();            
        end
        
        function computeFunctionAndGradient(obj)
            obj.computePhysicalData();
            obj.computeCCstar();
            
            
            obj.filterGradient();
            %obj.value = log10(obj.value);
            %obj.value = obj.value.^2;
            %obj.value = sum(obj.value);

        end
    end
end
