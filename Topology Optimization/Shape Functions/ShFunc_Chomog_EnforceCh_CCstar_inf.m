classdef ShFunc_Chomog_EnforceCh_CCstar_inf < ShFunc_Chomog_EnforceCh
    properties
        component
        epsilon
        initial_value = 1;
    end
    methods
        function obj=ShFunc_Chomog_EnforceCh_CCstar_inf(settings,n)
            obj.init(settings);
            obj.compute_Ch_star(settings.TOL, settings.selectiveC_Cstar);
            obj.component = n;
            obj.epsilon = settings.epsilon_isotropy_final;
        end
        function computeCostAndGradient(obj,x)
            obj.computePhysicalData(x);
            obj.computeCCstar(x);
            
            costval=obj.value;
            obj.value = (costval(obj.component))^2 - obj.epsilon;
            obj.gradient = 2*costval(obj.component)*obj.gradient(:,obj.component);
            obj.value = obj.value/obj.initial_value;
            obj.gradient = obj.gradient/obj.initial_value;
            obj.passFilter();
        end
    end
    
    methods(Access = public)
        function obj =  setEpsilon(obj,r)
            obj.epsilon = r;
        end
    end

end