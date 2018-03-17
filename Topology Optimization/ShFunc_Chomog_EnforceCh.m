classdef ShFunc_Chomog_EnforceCh< ShFunc_Chomog
     properties (Access = protected)  
         Ch_star
         selectiveC_Cstar
     end
    methods
        function obj=ShFunc_Chomog_EnforceCh(settings)     
            obj@ShFunc_Chomog(settings);
        end
        function obj = passFilter(obj)            
            mass=obj.filter.Msmooth;
            gradient=obj.filter.getP1fromP0(obj.gradient(:));
            gradient = mass*gradient;
            if isempty(obj.h_C_0)
                obj.h_C_0 = obj.value;
            end
%             obj.value = obj.value/abs(obj.h_C_0);
%             gradient=gradient/abs(obj.h_C_0);
%             obj.h_C_0 = costfunc;
                        
            obj.gradient = gradient;
        end
    end
end