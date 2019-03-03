classdef CantileverTriangleCoarse_Case_1_1_1 < TopOpt_Settings
    properties
 
    end
    methods
        function obj=CantileverTriangleCoarse_Case_1_1_1()
            obj.filter_settings.type='P1';
            obj.optimizer_settings=SLERP_Settings();
            obj.optimizer_settings.nconstr=1;
            obj.init_settings=CIRCLE_Settings();
%             obj.init_settings=DVar_Settings('circle');
%             obj.init_settings.case='circle';
            obj.init_settings.radius=1;
            
            
        end
    end
    
end