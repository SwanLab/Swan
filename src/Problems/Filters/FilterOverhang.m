classdef FilterOverhang < handle
    
    properties (Access = private)
        mesh
        trial
        k
        theta
        epsilon
        gamma
    end

    properties (Access = private)
        LHS
        prox
        rhsProx
        chiN
        chiNOld
    end

    methods (Access = public)
        function obj = FilterOverhang(cParams)
            obj.init(cParams);
            obj.createLHS();
        end

        function xF = compute(obj,fun,q)
            if isempty(obj.chiNOld)
                xF = obj.computeInitialGuess(fun,q);
            else
                xF = copy(obj.trial);
            end
            obj.chiN = obj.createRHSShapeFunction(fun,q);
            if isempty(obj.chiNOld) || norm(obj.chiNOld-obj.chiN)/norm(obj.chiN)>=1e-6
                obj.chiNOld      = obj.chiN;
                obj.updateProximal(xF);
                obj.updateRHSProx();
                xF.setFValues(obj.LHS\(obj.chiN + obj.rhsProx));
                obj.trial = xF;
            end
        end

        function obj = updateEpsilon(obj,epsilon)
            if obj.hasEpsilonChanged(epsilon)
                h           = obj.mesh.computeMeanCellSize();
                obj.epsilon = epsilon;
                obj.gamma   = (1/(2*epsilon^2))*(epsilon^2/h^2 - 1);
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
        end

        function xF = computeInitialGuess(obj,fun,q)
            s.mesh  = obj.mesh;
            s.trial = obj.trial;
            filter  = FilterLump(s);
            xF      = filter.compute(fun,q);
        end

        function RHS = createRHSShapeFunction(obj,fun,quadType)
            f   = @(v) DP(fun,v);
            RHS = IntegrateRHS(f,obj.trial,obj.trial.mesh,'Domain',quadType);
        end

        function updateProximal(obj,rho)
            gRho   = Grad(rho);
            gRhoN  = sqrt(DP(Grad(rho),Grad(rho)));
            dirDer = DP(gRho,obj.k);
            constr = dirDer - gRhoN.*cos(obj.theta);
            coef   = 1/(obj.gamma*obj.epsilon^2+1);
            obj.prox = coef.*gRho.*(constr>=0) + gRho.*(constr<0);
        end

        function updateRHSProx(obj)
            f   = @(v) DP(obj.prox,Grad(v));
            obj.rhsProx = (1/obj.gamma).*IntegrateRHS(f,obj.trial,obj.trial.mesh,'Domain');
        end

        function createLHS(obj)
            f = @(v,u) v.*u;
            M = IntegrateLHS(f,obj.trial,obj.trial,obj.mesh,'Domain',2);
            %M = diag(sum(M,1));
            g = @(v,u) DP(Grad(v),Grad(u));
            K = IntegrateLHS(g,obj.trial,obj.trial,obj.mesh,'Domain',2);
            obj.LHS = M + (1/obj.gamma).*K;
        end

        function itHas = hasEpsilonChanged(obj,eps)
            var   = abs(eps - obj.epsilon)/eps;
            itHas = var > 1e-15;
        end
    end
end