classdef ComplianceFromConstiutiveTensor < handle

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
            strain      = SymGrad(u);
            stress      = DDP(C,strain);
            dCompliance = DDP(strain,stress);
            J           = Integrator.compute(dCompliance,obj.mesh,obj.quadrature.order);
        end

    end

    methods (Static, Access = private)

        function dj = computeGradient(dC,u)
            strain  = SymGrad(u);
            dStress = DDP(dC,strain);
            dj      = -DDP(strain, dStress);
        end

    end

end
