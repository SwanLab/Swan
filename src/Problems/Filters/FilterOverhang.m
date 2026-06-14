classdef FilterOverhang < handle
    
    properties (Access = private)
        mesh
        trial
        k
        theta
        epsilon
        gamma
        gammaMax
        tol
    end

    properties (Access = private)
        M
        K
        LHS
        prox
        rhsProx
        chiN
        chiNOld
    end

    methods (Access = public)
        function obj = FilterOverhang(cParams)
            obj.init(cParams);
            obj.createMass();
            obj.createStiffness();
            obj.updateLHS();
        end

        function xF = compute(obj,fun,q)
            xF = obj.computeInitialGuess(fun,q);
            obj.chiN = obj.createRHSShapeFunction(fun,q);
            if isempty(obj.chiNOld) || norm(obj.chiNOld-obj.chiN)/norm(obj.chiN)>=1e-6
                obj.chiNOld = obj.chiN;
                obj.updateProximal(xF);
                obj.updateRHSProx();
                value0 = obj.computeRefCost(fun,xF);
                mOld  = obj.computeCost(fun,xF)/value0;
                delta = 1;
                xFOld = copy(xF);
                while (delta>=1e-6)
                    xF.setFValues(full(obj.LHS\(obj.chiN + obj.rhsProx)));
                    dxF = xF - xFOld;
                    mNew = obj.computeCost(fun,xF)/value0;
                    if (mNew - mOld)/Norm(dxF,'L2') <= 1e-2
                        delta = Norm(dxF,'L2')/Norm(xF,'L2');
                        obj.gamma = min(obj.gamma*1.2,obj.gammaMax);
                        obj.updateLHS();
                        obj.updateProximal(xF);
                        obj.updateRHSProx();
                        xFOld.setFValues(xF.fValues);
                        mOld = mNew;
                    else
                        obj.gamma = obj.gamma/2;
                        obj.updateLHS();
                        obj.updateProximal(xFOld);
                        obj.updateRHSProx();
                    end
                end
                obj.trial = xF;
            end
        end

        function obj = updateEpsilon(obj,epsilon)
            if obj.hasEpsilonChanged(epsilon)
                h            = obj.mesh.computeMeanCellSize();
                obj.epsilon  = epsilon;
                obj.gammaMax = (1/(2*epsilon^2))*(epsilon^2/((3*h)^2) - 1);
                obj.gamma    = obj.gammaMax;
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
            obj.gamma   = 0.001;
            obj.tol     = 1e-6;
        end

        function xF = computeInitialGuess(obj,fun,q)
            if isempty(obj.chiNOld)
                s.mesh  = obj.mesh;
                s.trial = obj.trial;
                filter  = FilterLump(s);
                xF      = filter.compute(fun,q);
            else
                xF = copy(obj.trial);
            end
        end

        function RHS = createRHSShapeFunction(obj,fun,quadType)
            f   = @(v) DP(fun,v);
            RHS = IntegrateRHS(f,obj.trial,obj.trial.mesh,'Domain',quadType);
        end

        function updateProximal(obj,rho)
            [gRho,constr] = obj.computeIntermediateOperators(rho);
            coef   = 1/(obj.gamma*obj.epsilon^2+1);
            pCon   = DomainFunction.create(@(xV) double(constr.evaluate(xV)>=0),obj.mesh,1);
            nCon   = DomainFunction.create(@(xV) double(constr.evaluate(xV)<0),obj.mesh,1);
            obj.prox = coef.*gRho.*pCon + gRho.*nCon;
        end

        function [gRho,constr] = computeIntermediateOperators(obj,rho)
            gRho   = Grad(rho);
            gRhoN  = sqrt(DP(gRho,gRho));
            dirDer = DP(gRho,obj.k);
            constr = dirDer./(gRhoN + obj.tol) - cos(obj.theta);
        end

        function updateRHSProx(obj)
            f   = @(v) DP(obj.prox,Grad(v));
            obj.rhsProx = (1/obj.gamma).*IntegrateRHS(f,obj.trial,obj.trial.mesh,'Domain');
        end

        function J = computeMinimumSquaresTerm(obj,chi,rho)
            int1 = Integrator.compute(rho.*rho,obj.mesh,2);
            int2 = -2*Integrator.compute(chi.*rho,obj.mesh,2);
            int3 = Integrator.compute(chi.*chi,obj.mesh,2);
            J    = 0.5*(int1+int2+int3);
        end

        function J = computeRegularizationTerm(obj,rho)
            [gRho,constr] = obj.computeIntermediateOperators(rho);
            e             = obj.epsilon;
            Z             = ConstantFunction.create(0,obj.mesh);
            sgn           = constr./(abs(constr) + 0.01*obj.tol);
            maxF          = max(Z,sgn);
            f             = DP(gRho,gRho);
            int           = f.*(maxF.^2);
            J             = (e^2/2)*Integrator.compute(int,obj.mesh,3);
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
            obj.M = diag(sum(Mraw,1));
        end

        function createStiffness(obj)
            g = @(v,u) DP(Grad(v),Grad(u));
            obj.K = IntegrateLHS(g,obj.trial,obj.trial,obj.mesh,'Domain',2);
        end

        function updateLHS(obj)
            obj.LHS = obj.M + (1/obj.gamma).*obj.K;
        end

        function itHas = hasEpsilonChanged(obj,eps)
            var   = abs(eps - obj.epsilon)/eps;
            itHas = var > 1e-15;
        end
    end
end