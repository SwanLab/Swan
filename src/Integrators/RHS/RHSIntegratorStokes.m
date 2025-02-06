classdef RHSIntegratorStokes < RHSIntegrator

    properties (Access = private)
        velocityFun
        pressureFun
        forcesFormula
    end

    methods (Access = public)

        function obj = RHSIntegratorStokes(cParams)
            obj.init(cParams);
        end

        function rhs = integrate(obj)
            rhs = obj.computeRHS();
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.mesh          = cParams.mesh;
            obj.velocityFun   = cParams.velocityFun;
            obj.pressureFun   = cParams.pressureFun;
            obj.forcesFormula = cParams.forcesFormula;
        end

        function RHS = computeRHS(obj)
            Fext = obj.computeVolumetricFext();
            g = obj.computeVelocityDivergence();
            RHS = [Fext; zeros(obj.pressureFun.nDofs,1)];
        end

        function Fext = computeVolumetricFext(obj)
            a.type = 'ShapeFunction';
            a.mesh = obj.mesh;
            a.quadType = 3;
            rhsI       = RHSIntegrator.create(a);
            test = LagrangianFunction.create(obj.mesh, 2, 'P2');
            Fext = rhsI.compute(obj.forcesFormula,test);
        end

        function g = computeVelocityDivergence(obj)
            g = zeros(obj.pressureFun.nDofs,1);
        end

    end

end