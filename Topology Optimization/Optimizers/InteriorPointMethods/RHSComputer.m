classdef RHSComputer < handle

    properties (Access = public)
    end
    properties (Access = protected)
        nX
        nSlack
        m
        cost
        constraint
        lowerZ
        upperZ
        lambda
        baseVariables
        invDiagdL
        invDiagdU
        e
    end

    methods (Access = public, Static)
        function obj = create(cParams)
            f = RHSFactory();
            obj = f.create(cParams);
        end

    end
    methods (Access = protected)
        function init(obj,cParams)
            obj.nX = cParams.nX;
            obj.nSlack = cParams.nSlack;
            obj.m = cParams.m;
            obj.cost = cParams.cost;
            obj.lowerZ = cParams.lowerZ;
            obj.upperZ = cParams.upperZ;
            obj.constraint = cParams.constraint;
            obj.lambda = cParams.lambda;
            obj.baseVariables = cParams.baseVariables;
            obj.invDiagdL = cParams.invDiagdL;
            obj.invDiagdU = cParams.invDiagdU;
            obj.e = cParams.e;
        end
    end
end