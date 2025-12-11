classdef ComplianceFromConstitutiveTensorThermoElastic < handle

    properties (Access = private)
        quadrature
    end

    properties (Access = private)
        mesh
        stateProblem
        adjointProblem
    end

    methods (Access = public)
        function obj = ComplianceFromConstitutiveTensorThermoElastic(cParams)
            obj.init(cParams);
            obj.createQuadrature();
            obj.createAdjointProblem(cParams);
        end

        function [J,dJ] = computeFunctionAndGradient(obj,C,dC,kappa, dkappa)
            [u,T] = obj.computeStateVariable(C,kappa);
            [p] = obj.computeAdjointVariable(kappa,C,u);
            J  = obj.computeFunction(C,u);
            dJ = obj.computeGradient(dC,u,dkappa,T,p);
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

        function createAdjointProblem(obj,cParams)
            % define bcs + use structure from thermal Problems
            
            obj.adjointProblem = ThermalProblem(s); 
            
        end

        function [u,T] = computeStateVariable(obj,C,kappa)
            obj.stateProblem.updateMaterial(C,kappa);
            obj.stateProblem.solve();
            u = obj.stateProblem.uFun;
            T = obj.stateProblem.tFun;
        end

        function [p] = computeAdjointVariable(kappa,C,u)
            I = ConstantFunction.create(eye(2),obj.mesh);
            beta = obj.alpha.*DDP(C,I);    
            newSource =  DDP(beta,SymGrad(u));        % updating RHS
            obj.adjointProblem.updateSource(newSource);
            obj.adjointProblem.solve(kappa);
            p = obj.adjointProblem.pFun;
        end

        function J = computeFunction(obj,C,u)
            dCompliance = ElasticEnergyDensity(C,u);
            J           = Integrator.compute(dCompliance,obj.mesh,obj.quadrature.order);
        end
    end

    methods (Static, Access = private)
        function dj = computeGradient(dC,u,dkappa,T)
            nDesVar = length(dC);
            dj      = cell(nDesVar,1);
            for i = 1:nDesVar
                strain  = SymGrad(u);
                dStress = DDP(dC{i},strain);
                I = ConstantFunction.create(eye(2),obj.mesh);
                dbeta = obj.alpha.*DDP(dC{i},I);       
                dj{i}   = -0.5.*DDP(strain, dStress) + DDP(dbeta.*Grad(T),SymGrad(u)) + DP(Grad(T),dkappa.*Grad(p)); %T*dkappa12*u + T*dkappa*p + (-dQ*p+df*u)
            end
        end
    end
end