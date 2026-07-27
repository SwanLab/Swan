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
        rieszEps
    end

    properties (Access = private)
        M
        K
        LHS
        LHSRiesz
        s
        chiN
        chiNOld
        rhsProx
        rhsDualNorm
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
            if isempty(obj.chiNOld) || norm(obj.chiNOld-obj.chiN)/norm(obj.chiN)>=1e-6
                xF = obj.computeInitialGuess(fun,q);
                obj.chiNOld = obj.chiN;
                obj.updateProximalGrad(xF);
                obj.updateRHSProx(xF,q);
                obj.updateRHSDualNorm();
                dJ0 = obj.computeCostGradient(xF);
                error0 = Norm(dJ0,'H1',obj.rieszEps);
                error = inf;
                xFOld = copy(xF);
                iter = 1;
                while (error>=obj.tol && error0>=obj.tol && iter<=10000)
                    t  = obj.tau;
                    dx = full(obj.LHS\(t.*(obj.chiN - obj.rhsProx - obj.rhsDualNorm)));
                    xF.setFValues(xFOld.fValues + dx);
                    obj.updateProximalGrad(xF);
                    obj.updateRHSDualNorm();
                    dJ = obj.computeCostGradient(xF);
                    error = Norm(dJ,'H1',obj.rieszEps);
                    if error < error0 || obj.tau<1e-5
                        obj.tau = min(obj.tau*1.1,obj.tauMax);
                        obj.updateLHS();
                        obj.updateRHSProx(xF,q);
                        xFOld.setFValues(xF.fValues);
                        error0 = error;
                        iter = iter + 1;
                    else
                        obj.tau = obj.tau/1.5;
                        obj.updateLHS();
                        obj.updateProximalGrad(xFOld);
                        obj.updateRHSProx(xFOld,q);
                        obj.updateRHSDualNorm();
                    end
%                     if obj.tau<=1.2e-5 && abs(error-error0)<obj.tol
%                         break;
%                     end
                end
                obj.trial = xF;
            end
            xF = obj.trial;
        end

        function obj = updateEpsilon(obj,epsilon)
            if obj.hasEpsilonChanged(epsilon)
                obj.epsilon  = epsilon;
                obj.tauMax   = 10;
                obj.tau      = 0.1;
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
            obj.tauMax    = 10;
            obj.tau       = 0.1;
            obj.tolPrimal = 1e-12;
            obj.tol       = cParams.tol;
            obj.theta0    = deg2rad(2*15);
            obj.rieszEps  = 8*obj.mesh.computeMeanCellSize();
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

        function updateProximalGrad(obj,rho)
            obj.s = Grad(rho);
        end

        function phi = computeAngleRespectDiamondAxis(obj,gRho)
            gRhoN   = sqrt(DP(gRho,gRho));
            dirDer  = DP(gRho,obj.k);
            cosPhi  = dirDer./(gRhoN + obj.tolPrimal);
            phi     = acos(cosPhi);
        end

        function updateRHSProx(obj,xF,q)
            obj.rhsProx = obj.createRHSShapeFunction(xF,q);
        end

        function updateRHSDualNorm(obj)
            phi = obj.computeAngleRespectDiamondAxis(obj.s);
            th0 = obj.theta0;
            h   = obj.mesh.computeMeanCellSize();
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
            obj.rhsDualNorm = (e^2/2).*IntegrateRHS(der,obj.trial,obj.mesh,'Domain',4);
        end

        function dJ = computeCostGradient(obj,xF)
            rhsRhoe = obj.M*xF.fValues;
            rhs = rhsRhoe - obj.chiN + obj.rhsDualNorm;
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
            t = obj.tau;
            e = obj.rieszEps;
            obj.LHS = (1+t).*obj.M + e^2.*obj.K;
        end

        function createLHSRiesz(obj)
            e = obj.rieszEps;
            obj.LHSRiesz = decomposition(obj.M + e^2*obj.K);
        end

        function itHas = hasEpsilonChanged(obj,eps)
            var   = abs(eps - obj.epsilon)/eps;
            itHas = var > 1e-15;
        end
    end
end