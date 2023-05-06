classdef SymmetricLHSComputer < LHSComputer

    properties (Access = public)
        LHS
    end
    properties (Access = private)
    end

    methods (Access = public)
        function obj = SymmetricLHSComputer(cParams)
            obj.init(cParams);
        end
        function compute(obj)
            obj.computeLHS();
        end
    end

    methods (Access = private)
        function computeLHS(obj)
            obj.LHS = [obj.H,obj.constraint.gradient;obj.constraint.gradient',zeros(obj.m,obj.m)];
        end
    end
end