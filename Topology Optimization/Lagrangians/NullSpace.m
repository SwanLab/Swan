classdef NullSpace < ObjectiveFunction

    properties (Access = private)
        constraintModifier
    end

    

    methods (Access = public)

        function obj = NullSpace(cParams)
            obj.init(cParams);
            obj.createConstraintModifier(cParams);
        end

        function updateBecauseOfPrimal(obj)
            obj.cost.computeFunctionAndGradient();
            obj.constraint.computeFunctionAndGradient();
            obj.computeFunction();
            obj.computeGradient();
        end

        function obj = updateBecauseOfDual(obj)
            obj.modifyInactiveConstraints();
            obj.computeFunction();
            obj.computeGradient();
        end

        function createConstraintModifier(obj,cParams)
            s.type         = cParams.constraintCase;
            s.dualVariable = cParams.dualVariable;
            s.constraint   = cParams.constraint;
            icm = InactiveConstraintsModifier.create(s);
            obj.constraintModifier = icm;
        end

        function computeFunction(obj)
            l         = obj.dualVariable.value;
            C         = obj.constraint.value;
            DC        = obj.constraint.gradient';
            J         = obj.cost.value;
            alphaJ    = 0.001;
            alphaC    = 0.001;
            S         = (DC*DC')^-1;
            v         = alphaJ*(J + l*C) + alphaC/2*C'*S*C;
            obj.value = v;
        end

        function computeGradient(obj)
            C         = obj.constraint.value;
            Jgrad     = obj.cost.gradient;
            alphaJ    = 0.001;
            alphaC    = 0.001;
            DC        = obj.constraint.gradient';
            S         = (DC*DC')^-1;
            A         = DC'*S*DC;
            I         = ones(size(A,1),1);
            epsilonJ  = (I - A)*Jgrad;
            epsilonC  = DC'*S*C;
            g         = alphaJ*epsilonJ + alphaC*epsilonC;
            obj.gradient = g;
        end

    end
    
end