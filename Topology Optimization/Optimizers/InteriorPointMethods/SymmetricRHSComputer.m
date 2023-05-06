classdef SymmetricRHSComputer < RHSComputer

    properties (Access = public)
        RHS
    end
    properties (Access = private)
    end

    methods (Access = public)
        function obj = SymmetricRHSComputer(cParams)
            obj.init(cParams)
        end
        function compute(obj)
            obj.computeRHS();
        end
    end
    methods (Access = private)
        function computeRHS(obj)
            obj.RHS(1:obj.nX + obj.nSlack,1) = obj.cost.gradient'...
                - obj.lowerZ' + obj.upperZ' + obj.constraint.gradient*obj.lambda' +...
                obj.lowerZ'-obj.baseVariables.mu*obj.invDiagdL*obj.e -obj.upperZ'...
                +obj.baseVariables.mu*obj.invDiagdU*obj.e;
            obj.RHS(obj.nX + obj.nSlack + 1:obj.nX + obj.nSlack + obj.m,1)...
                = obj.constraint.value(1:obj.m)';
        end
    end

end