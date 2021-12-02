classdef ShFunc_Chomog_EnforceCh_CCstar_eq < ShFunc_Chomog_EnforceCh
    properties
        component
        initial_value = 1;
    end
    methods
        function obj=ShFunc_Chomog_EnforceCh_CCstar_eq(settings,n)
            obj.init(settings);
            obj.compute_Ch_star(settings.TOL, settings.selectiveC_Cstar);
            obj.component = n;
        end
        
        function computeCostAndGradient(obj,x)
            obj.computePhysicalData(x);
            obj.computeCCstar(x);
            
            obj.value = obj.value(obj.component);
            obj.gradient = obj.gradient(:,obj.component);
            obj.value = obj.value/obj.initial_value;
            obj.gradient = obj.gradient/obj.initial_value;
            obj.passFilter();
        end
    end
end

