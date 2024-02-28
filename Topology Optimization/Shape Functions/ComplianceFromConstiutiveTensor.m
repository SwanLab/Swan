classdef ComplianceFromConstiutiveTensor < handle

    properties (Access = public)

    end

    properties (Access = private)
        quadrature
    end

    properties (Access = private)
        mesh
        stateProblem
    end

    methods (Access = public)

        function obj = ComplianceFromConstiutiveTensor(cParams)
            obj.init(cParams);
            obj.createQuadrature();
        end

        function [J,dJ] = computeFunctionAndGradient(obj,C,dC)
            u  = obj.computeStateVariable(C);
            J  = obj.computeFunction(C,u);
            dJ = obj.computeGradient(dC,u);
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.mesh         = cParams.mesh;
            obj.stateProblem = cParams.stateProblem;
        end

        function createQuadrature(obj)
            quad = Quadrature.set(obj.mesh.type);
            quad.computeQuadrature('QUADRATIC');
            obj.quadrature = quad;
        end

        function u = computeStateVariable(obj,C)
            obj.stateProblem.updateMaterial(C);
            obj.stateProblem.solve();
            u = obj.stateProblem.uFun;
        end

        function J = computeFunction(obj,C,u)
            strain = SymGrad(u);
            stress = DDP(C,strain);
            int    = Integrator.create('ScalarProduct',obj.mesh,obj.quadrature.order);
            J      = int.compute(strain,stress);
        end
  
        function dj = computeGradient(obj,C,u)
            eu2 = SymGrad(u);
            % dStr = DDP(dC, eu2);
            dj = -DDP(eu2, C, eu2); % !!
        end

        function fd = createGaussFunction(obj,f)
            m  = obj.mesh;
            q  = obj.quadrature;
            fd = FGaussDiscontinuousFunction.create(f,m,q);
        end

    end

end
