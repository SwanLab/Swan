classdef Projector < handle

    properties (Access = protected)

    end

    methods (Static, Access = public)

        function obj = create(cParams)
            p = ProjectorFactory();
            obj = p.create(cParams);
        end

    end

    methods (Access = protected)

    end

end