classdef Monitoring < handle

    methods (Access = public, Static)
        function obj = create(cParams)
            m = MonitoringFactory();
            obj = m.create(cParams);
        end
    end
end