classdef FilterRegularizedDiamond < handle
    
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
        theta0
    end

    properties (Access = private)
        M
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
        function obj = FilterRegularizedDiamond(cParams)
            obj.init(cParams);
            obj.createMass();
            obj.createStiffness();
            obj.updateLHS();
            obj.createLHSRiesz();
        end

        function xF = compute(obj,fun,q)
            obj.chiN = obj.createRHSShapeFunction(fun,q);
            h = obj.mesh.computeMeanCellSize();
            if isempty(obj.chiNOld) || norm(obj.chiNOld-obj.chiN)/norm(obj.chiN)>=1e-6
                xF = obj.computeInitialGuess(fun,q);
                obj.chiNOld = obj.chiN;
                obj.updateProximal(xF);
                obj.updateRHSProx();
                obj.updateRHSDualNorm();
                dJ0 = obj.computeCostGradient(xF);
                error0 = Norm(dJ0,'H1',8*h);
                error = inf;
                xFOld = copy(xF);
                iter = 1;
                while (error>=obj.tol && error0>=obj.tol && iter<=10000)
                    xF.setFValues(full(obj.LHS\(obj.chiN + obj.rhsProx + obj.rhsDualNorm)));
                    obj.updateProximal(xF);
                    obj.updateRHSDualNorm();
                    dJ = obj.computeCostGradient(xF);
                    error = Norm(dJ,'H1',8*h);
                    if error <= (1.000001)*error0 %&& abs(Norm(xF-xFOld,'L2'))/abs(Norm(ConstantFunction.create(1,obj.mesh),'L2')) <= 0.2
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
                obj.tau      = 100;
                obj.updateLHS();
                obj.chiNOld = [];
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
            obj.tauMax    = 10000;
            obj.tau       = 100;
            obj.tolPrimal = 1e-12;
            obj.tol       = cParams.tol;
            obj.theta0    = deg2rad(2*15);
        end

        function xF = computeInitialGuess(obj,fun,q)
            if sum(obj.trial.fValues) == 0
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

        function phi = computeAngleRespectDiamondAxis(obj,gRho)
            gRhoN   = sqrt(DP(gRho,gRho));
            dirDer  = DP(gRho,obj.k);
            cosPhi  = dirDer./(gRhoN + obj.tolPrimal);
            phi     = acos(cosPhi);
        end

        function updateRHSProx(obj)
            f   = @(v) DP(obj.s,Grad(v));
            obj.rhsProx = (1/obj.tau).*IntegrateRHS(f,obj.trial,obj.trial.mesh,'Domain');
        end

        function updateRHSDualNorm(obj)
            phi = obj.computeAngleRespectDiamondAxis(obj.s);
            th0 = obj.theta0;
            h   = 1.5*obj.mesh.computeMeanCellSize();
            e   = obj.epsilon;

            heav = tanh((phi - obj.theta)./(th0/4));
            coef = 1+h/e - (1-h/e).*heav;

            der1 = @(u) ((1 + h/e)/2).*DP(obj.s,Grad(u));
            der2 = @(u) ((1 - h/e)/2).*DP(obj.s.*heav,Grad(u));

            l2N  = sqrt(DP(obj.s,obj.s)) + obj.tolPrimal;
            cosPhiDer = @(u) DP(Grad(u),obj.k)./l2N - DP((obj.s./l2N).*DP(obj.s./l2N,Grad(u))./l2N,obj.k);
            phiDer = @(u) - cosPhiDer(u)./(abs(sin(phi)) + obj.tolPrimal);
            der3 = @(u) ((1 - h/e)/2)*DP(obj.s,obj.s).*((sech((phi - obj.theta)./(th0/4))).^2).*phiDer(u)./(th0/4);

            der = @(u) coef.*(der1(u) - der2(u) - der3(u));
            obj.rhsDualNorm = (-e^2/2).*IntegrateRHS(der,obj.trial,obj.mesh,'Domain',4);
        end

        function dJ = computeCostGradient(obj,xF)
            rhsRhoe = obj.M*xF.fValues;
            rhs = rhsRhoe - obj.chiN - obj.rhsDualNorm;
            fVal = obj.LHSRiesz\rhs;
            dJ = copy(obj.trial);
            dJ.setFValues(fVal);
        end

        function createMass(obj)
            f = @(v,u) v.*u;
            obj.M = IntegrateLHS(f,obj.trial,obj.trial,obj.mesh,'Domain',2);
        end

        function createStiffness(obj)
            g = @(v,u) DP(Grad(v),Grad(u));
            obj.K = IntegrateLHS(g,obj.trial,obj.trial,obj.mesh,'Domain',2);
        end

        function updateLHS(obj)
            obj.LHS = obj.M + (1/obj.tau).*obj.K;
        end

        function createLHSRiesz(obj)
            h = obj.mesh.computeMeanCellSize();
            obj.LHSRiesz = decomposition(obj.M + (8*h)^2*obj.K);
        end

        function itHas = hasEpsilonChanged(obj,eps)
            var   = abs(eps - obj.epsilon)/eps;
            itHas = var > 1e-15;
        end
    end
end