classdef IPMSymmetricRHSComputer < handle

    properties (Access = private)
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

    properties (Access = private)
        diagonaldL
        diagonaldU
    end

    methods (Access = public)
        function obj = IPMSymmetricRHSComputer(cParams)
            obj.init(cParams)
        end
        function RHS = compute(obj)
            RHS= obj.computeRHS();
        end
    end
    methods (Access = private)

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

        function RHS = computeRHS(obj)
            RHS(1:obj.nX + obj.nSlack,1) = obj.cost.gradient'...
                - obj.lowerZ' + obj.upperZ' + obj.constraint.gradient*obj.lambda' +...
                obj.lowerZ'-obj.baseVariables.mu.*(obj.invDiagdL*obj.e) -obj.upperZ'...
                +obj.baseVariables.mu.*(obj.invDiagdU*obj.e);
            RHS(obj.nX + obj.nSlack + 1:obj.nX + obj.nSlack + obj.m,1)...
                 = obj.constraint.value(1:obj.m)';
        end
    end

end