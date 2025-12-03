classdef ComplianceFunctionalRadius < handle

    properties (Access = private)
        quadrature
    end

    properties (Access = private)
        mesh
        stateProblem
    end

    methods (Access = public)
        function obj = ComplianceFunctionalRadius(cParams)
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
            quad = Quadrature.create(obj.mesh,2);
            obj.quadrature = quad;
        end

        function u = computeStateVariable(obj,C)
            obj.stateProblem.updateMaterial(C);
            obj.stateProblem.solve();
            u = obj.stateProblem.uFun;
        end

        function J = computeFunction(obj,C,u)
            dCompliance = ElasticEnergyDensity(C,u);
            J           = Integrator.compute(dCompliance,obj.mesh,obj.quadrature.order);
        end
    end

    methods (Static, Access = private)
        function dj = computeGradient(dC,u)
            nDesVar = length(dC);
            dj      = cell(nDesVar,1);
            for i = 1:nDesVar
                strain  = SymGrad(u);
                dStress = DDP(dC{i},strain);
                dj{i}   = -0.5.*DDP(strain, dStress);
            end
        end
    end
end