classdef pyAMG < handle

    methods (Static, Access = public)
        function obj = create(s)
            switch s.type
                case 'ELASTIC'
                    obj = SmoothedAggregation(s);
            end
        end
    end
end