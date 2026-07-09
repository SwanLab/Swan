classdef FilterDiamond < handle
    
    properties (Access = private)
        mesh
        trial
        k
        theta
        epsilon
        tau
        tauMax
        tol
    end

    properties (Access = private)
        M
        K
        LHS
        s
        rhsProx
        rhsDualNorm
        chiN
        chiNOld
    end

    methods (Access = public)
        function obj = FilterDiamond(cParams)
            obj.init(cParams);
            obj.createMass();
            obj.createStiffness();
            obj.updateLHS();
        end

        function xF = compute(obj,fun,q)
            obj.chiN = obj.createRHSShapeFunction(fun,q);
            if isempty(obj.chiNOld) || norm(obj.chiNOld-obj.chiN)/norm(obj.chiN)>=1e-6
                xF = obj.computeInitialGuess(fun,q);
                obj.chiNOld = obj.chiN;
                obj.updateProximal(xF);
                obj.updateRHSProx();
                obj.updateRHSDualNorm();
                value0 = obj.computeRefCost(fun,xF);
                mOld  = obj.computeCost(fun,xF)/value0;
                delta = 1;
                xFOld = copy(xF);
                iter = 1;
                while (delta>=1e-6 && iter<=10000)
                    xF.setFValues(full(obj.LHS\(obj.chiN + obj.rhsProx + obj.rhsDualNorm)));
                    mNew = obj.computeCost(fun,xF)/value0;
                    if mNew - mOld <= 1e-10
                        delta = abs(mNew - mOld);
                        obj.tau = min(obj.tau*1.2,obj.tauMax);
                        obj.updateLHS();
                        obj.updateProximal(xF);
                        obj.updateRHSProx();
                        obj.updateRHSDualNorm();
                        xFOld.setFValues(xF.fValues);
                        mOld = mNew;
                    else
                        obj.tau = obj.tau/2;
                        obj.updateLHS();
                        obj.updateProximal(xFOld);
                        obj.updateRHSProx();
                        obj.updateRHSDualNorm();
                    end
                    iter = iter + 1;
                end
                obj.trial = xF;
            end
            xF = obj.trial;
        end

        function obj = updateEpsilon(obj,epsilon)
            if obj.hasEpsilonChanged(epsilon)
                h            = obj.mesh.computeMeanCellSize();
                obj.epsilon  = epsilon;
                obj.tauMax   = 10000;
                obj.tau      = 1;
                obj.updateLHS();
            end
        end
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.trial   = LagrangianFunction.create(cParams.mesh, cParams.trial.ndimf, cParams.trial.order);
            obj.mesh    = cParams.mesh;
            obj.k       = cParams.senseVector;
            obj.theta   = deg2rad(90 - cParams.ovAngleDeg);
            obj.epsilon = cParams.mesh.computeMeanCellSize();
            obj.tau     = 1;
            obj.tol     = 1e-6;
        end

        function xF = computeInitialGuess(obj,fun,q)
            if isempty(obj.chiNOld)
                sF.mesh = obj.mesh;
                sF.trial = obj.trial;
                sF.filterType = 'PDE';
                filter = Filter.create(sF);
                filter.updateEpsilon(3*obj.mesh.computeMeanCellSize());
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
            cosPhi  = dirDer./(gRhoN + obj.tol);
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

            l2N  = sqrt(DP(obj.s,obj.s)) + obj.tol;
            cosPhiDer = @(u) DP(Grad(u),obj.k)./l2N - DP((obj.s./l2N).*DP(obj.s./l2N,Grad(u))./l2N,obj.k);
            sinPhiDer = @(u) -cos(phi).*cosPhiDer(u)./(sqrt(1-(cos(phi)).^2) + obj.tol);
            cosPhiThetaDer = @(u) cosPhiDer(u).*cos(obj.theta) + sinPhiDer(u).*sin(obj.theta);

            Z    = ConstantFunction.create(0,obj.mesh);
            maxF = max(Z,cos(phi-obj.theta));
            derOut = @(u) DP(2.*obj.s.*maxF.^2,Grad(u)) + DP(2.*obj.s.*maxF,obj.s.*cosPhiThetaDer(u));

            der = @(u) 2.*DP(obj.s,Grad(u)).*inDiam + derOut(u).*outDiam;
            e = obj.epsilon;
            obj.rhsDualNorm = (-e^2/2).*IntegrateRHS(der,obj.trial,obj.mesh,'Domain',3);
        end

        function J = computeMinimumSquaresTerm(obj,chi,rho)
            int1 = Integrator.compute(rho.*rho,obj.mesh,2);
            int2 = -2*Integrator.compute(chi.*rho,obj.mesh,2);
            int3 = Integrator.compute(chi.*chi,obj.mesh,2);
            J    = 0.5*(int1+int2+int3);
        end

        function J = computeRegularizationTerm(obj,rho)
            [phi,inDiam,outDiam] = obj.computeIntermediateOperators(Grad(rho));
            l2N  = DP(Grad(rho),Grad(rho));
            Z    = ConstantFunction.create(0,obj.mesh);
            maxF = max(Z,cos(phi-obj.theta)).^2;
            int  = l2N.*inDiam + l2N.*maxF.*outDiam;
            e    = obj.epsilon;
            J    = (e^2/2)*Integrator.compute(int,obj.mesh,3);
        end

        function J = computeRefCost(obj,chi,rho)
            int1 = Integrator.compute(rho.*rho,obj.mesh,2);
            int2 = 2*Integrator.compute(chi.*rho,obj.mesh,2);
            int3 = Integrator.compute(chi.*chi,obj.mesh,2);
            J    = 0.5*(int1+int2+int3);
        end

        function J = computeCost(obj,fun,rho)
            J1 = obj.computeMinimumSquaresTerm(fun,rho);
            J2 = obj.computeRegularizationTerm(rho);
            J  = J1 + J2;
        end

        function createMass(obj)
            f     = @(v,u) v.*u;
            Mraw  = IntegrateLHS(f,obj.trial,obj.trial,obj.mesh,'Domain',2);
            obj.M = Mraw;
        end

        function createStiffness(obj)
            g = @(v,u) DP(Grad(v),Grad(u));
            obj.K = IntegrateLHS(g,obj.trial,obj.trial,obj.mesh,'Domain',2);
        end

        function updateLHS(obj)
            obj.LHS = obj.M + (1/obj.tau).*obj.K;
        end

        function itHas = hasEpsilonChanged(obj,eps)
            var   = abs(eps - obj.epsilon)/eps;
            itHas = var > 1e-15;
        end
    end
end