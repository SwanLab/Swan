classdef ShFunc_Chomog_EnforceCh_CCstar_L2 < ShFunc_Chomog_EnforceCh
    properties
        initial_value
        epsilon
    end
    methods
        function obj=ShFunc_Chomog_EnforceCh_CCstar_L2(settings)
            obj@ShFunc_Chomog_EnforceCh(settings);
            obj.compute_Ch_star(settings.TOL, settings.selectiveC_Cstar);
            obj.epsilon = settings.epsilon_isotropy_final;
        end
        
        function computef(obj,x)
            obj.computePhysicalData(x);
            obj.computeCCstar(x);
            ct=1;
            
            %Gradient
            for i=1:6
                obj.gradient(:,i) = 2*obj.value(i)*obj.gradient(:,i);
            end
            obj.passFilter;
            obj.gradient = sum(obj.gradient,2);
            
            %Cost
            obj.value = obj.value.^2;
            if isempty(obj.initial_value)
                obj.initial_value = ct*sum(obj.value);
            end
            
            obj.value = ct*sum(obj.value)/obj.initial_value - obj.epsilon;            
            obj.gradient=ct*obj.gradient/obj.initial_value;
        end
    end
end
