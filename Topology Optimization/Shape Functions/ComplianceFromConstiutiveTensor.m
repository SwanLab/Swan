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

        function [J,dJ] = computeFunctionAndGradient(obj,C)
            u  = obj.computeStateVariable(C);
            J  = obj.computeFunction(C,u);
            dJ = obj.computeGradient(C,u);
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
  
        function g = computeGradient(obj,dC,u)
            g = obj.computeDJ(dC,u);
        end

        function dj = computeDJ(obj,dC,u)
            xV = obj.quadrature.posgp;
            eu2 = SymGrad(u);
            ngaus = size(xV,2);
            nelem = obj.mesh.nelem;
            % dStr = DDP(dC, eu2);
            dj = -DDP(eu2, dC, eu2); % !!
            dj = squeezeParticular(dj, 1, [1 1 ngaus nelem]);
        end

        function fd = createGaussFunction(obj,f)
            m  = obj.mesh;
            q  = obj.quadrature;
            fd = FGaussDiscontinuousFunction.create(f,m,q);
        end

    end

end
