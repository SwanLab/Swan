classdef FilterParabolicDiamond < handle
    
    properties (Access = private)
        mesh
        trial
        k
        theta
        epsilon
        tau
        tauMax
        tolPrimal
        tol
        parCoef
    end

    properties (Access = private)
        M
        Mlump
        K
        LHS
        LHSRiesz
        s
        rhsProx
        rhsDualNorm
        chiN
        chiNOld
    end

    methods (Access = public)
        function obj = FilterParabolicDiamond(cParams)
            obj.init(cParams);
            obj.createMass();
            obj.createStiffness();
            obj.updateLHS();
            obj.createLHSRiesz();
        end

        function xF = compute(obj,fun,q)
            obj.chiN = obj.createRHSShapeFunction(fun,q);
            if isempty(obj.chiNOld) || norm(obj.chiNOld-obj.chiN)/norm(obj.chiN)>=1e-6
                xF = obj.computeInitialGuess(fun,q);
                obj.chiNOld = obj.chiN;
                obj.updateProximal(xF);
                obj.updateRHSProx();
                obj.updateRHSDualNorm();
                dJ0 = obj.computeCostGradient(xF);
                error0 = Norm(dJ0,'L2');
                error = inf;
                xFOld = copy(xF);
                iter = 1;
                while (error>=obj.tol && error0>=obj.tol && iter<=10000)
                    xF.setFValues(full(obj.LHS\(obj.chiN + obj.rhsProx + obj.rhsDualNorm)));
                    obj.updateProximal(xF);
                    obj.updateRHSDualNorm();
                    dJ = obj.computeCostGradient(xF);
                    error = Norm(dJ,'L2');
                    if error <= 1.05*error0
                        obj.tau = min(obj.tau*1.2,obj.tauMax);
                        obj.updateLHS();
                        obj.updateRHSProx();
                        xFOld.setFValues(xF.fValues);
                        error0 = error;
                        iter = iter + 1;
                    else
                        obj.tau = obj.tau/2;
                        obj.updateLHS();
                        obj.updateProximal(xFOld);
                        obj.updateRHSProx();
                        obj.updateRHSDualNorm();
                    end
                end
                obj.trial = xF;
            end
            xF = obj.trial;
        end

        function obj = updateEpsilon(obj,epsilon)
            if obj.hasEpsilonChanged(epsilon)
                obj.epsilon  = epsilon;
                obj.tauMax   = 10000;
                obj.tau      = 1;
                obj.updateLHS();
            end
        end
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.trial     = LagrangianFunction.create(cParams.mesh, cParams.trial.ndimf, cParams.trial.order);
            obj.mesh      = cParams.mesh;
            obj.k         = cParams.senseVector;
            obj.theta     = deg2rad(90 - cParams.ovAngleDeg);
            obj.epsilon   = cParams.mesh.computeMeanCellSize();
            obj.tau       = 1;
            obj.tolPrimal = 1e-6;
            obj.tol       = cParams.tol;
            obj.parCoef   = (180/(20*pi))^2;
        end

        function xF = computeInitialGuess(obj,fun,q)
            if isempty(obj.chiNOld)
                sF.mesh = obj.mesh;
                sF.trial = obj.trial;
                sF.filterType = 'PDE';
                filter = Filter.create(sF);
                filter.updateEpsilon(obj.mesh.computeMeanCellSize());
                xF = filter.compute(fun,q);
            else
                xF = copy(obj.trial);
            end
        end

        function RHS = createRHSShapeFunction(obj,fun,quadType)
            f   = @(v) DP(fun,v);
            RHS = IntegrateRHS(f,obj.trial,obj.trial.mesh,'Domain',quadType);
        end

        function updateProximal(obj,rho)
            obj.s = Grad(rho);
        end

        function [phi,inDiam,outDiam] = computeIntermediateOperators(obj,gRho)
            gRhoN   = sqrt(DP(gRho,gRho));
            dirDer  = DP(gRho,obj.k);
            cosPhi  = dirDer./(gRhoN + obj.tolPrimal);
            phi     = acos(cosPhi);
            inDiam  = DomainFunction.create(@(xV) double(phi.evaluate(xV)<=obj.theta),obj.mesh,1);
            outDiam = DomainFunction.create(@(xV) double(phi.evaluate(xV)>obj.theta),obj.mesh,1);
        end

        function updateRHSProx(obj)
            f   = @(v) DP(obj.s,Grad(v));
            obj.rhsProx = (1/obj.tau).*IntegrateRHS(f,obj.trial,obj.trial.mesh,'Domain');
        end

        function updateRHSDualNorm(obj)
            [phi,inDiam,outDiam] = obj.computeIntermediateOperators(obj.s);

            coef = obj.parCoef;

            l2N  = sqrt(DP(obj.s,obj.s)) + obj.tolPrimal;
            cosPhiDer = @(u) DP(Grad(u),obj.k)./l2N - DP((obj.s./l2N).*DP(obj.s./l2N,Grad(u))./l2N,obj.k);

            Z    = ConstantFunction.create(0,obj.mesh);
            maxF = max(Z,1-coef.*(phi-obj.theta).^2);
            derOut = @(u) DP(2.*obj.s.*maxF.^2,Grad(u)) + DP(2.*obj.s.*maxF,2.*obj.s.*coef.*(phi-obj.theta).*cosPhiDer(u)./abs(sin(phi)));

            der = @(u) 2.*DP(obj.s,Grad(u)).*inDiam + derOut(u).*outDiam;
            e = obj.epsilon;
            obj.rhsDualNorm = (-e^2/2).*IntegrateRHS(der,obj.trial,obj.mesh,'Domain',5);
        end

        function dJ = computeCostGradient(obj,xF)
            rhsRhoe = obj.Mlump*xF.fValues;
            rhs = rhsRhoe - obj.chiN - obj.rhsDualNorm;
            fVal = obj.LHSRiesz\rhs;
            dJ = copy(obj.trial);
            dJ.setFValues(fVal);
        end

        function createMass(obj)
            f = @(v,u) v.*u;
            obj.M = IntegrateLHS(f,obj.trial,obj.trial,obj.mesh,'Domain',2);
            obj.Mlump = diag(sum(obj.M,1));
        end

        function createStiffness(obj)
            g = @(v,u) DP(Grad(v),Grad(u));
            obj.K = IntegrateLHS(g,obj.trial,obj.trial,obj.mesh,'Domain',2);
        end

        function updateLHS(obj)
            obj.LHS = obj.Mlump + (1/obj.tau).*obj.K;
        end

        function createLHSRiesz(obj)
            h = obj.mesh.computeMeanCellSize();
            obj.LHSRiesz = decomposition(obj.M + h^2*obj.K);
        end

        function itHas = hasEpsilonChanged(obj,eps)
            var   = abs(eps - obj.epsilon)/eps;
            itHas = var > 1e-15;
        end
    end
end