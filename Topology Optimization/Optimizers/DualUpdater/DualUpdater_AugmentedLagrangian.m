classdef DualUpdater_AugmentedLagrangian < DualUpdater
    
    properties (Access = private)
       augmentedLagrangian
       constraint
    end
        
    methods (Access = public)
        
        function obj = DualUpdater_AugmentedLagrangian(cParams)
            obj.init(cParams)
            obj.constraint = cParams.constraint;
            obj.augmentedLagrangian = cParams.augmentedLagrangian;
        end
        
        function updateDualVariable(obj)
            l   = obj.dualVariable.value;
            rho = obj.augmentedLagrangian.penalty;
            c   = obj.constraint.value';
            l = l + rho.*c;
            obj.dualVariable.value = l;
        end
        
    end
    
    methods (Access = private)
           
    end
    
  
end