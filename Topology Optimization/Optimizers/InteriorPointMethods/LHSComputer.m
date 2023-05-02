classdef LHSComputer < handle

    properties (Access = public)
    end
    properties (Access = protected)
        funcType
        H
        constraint
        m
    end

    methods (Access = public, Static)
        function obj = create(cParams)
            f = LHSFactory();
            obj = f.create(cParams);
        end
    end
    methods (Access = protected)
        function init(obj,cParams)
            obj.H = cParams.H;
            obj.constraint = cParams.constraint;
            obj.m = cParams.m;
        end
    end

end