classdef LHSIntegratorStokes < handle

    properties (Access = private)
        dt
        mesh
        velocityFun
        pressureFun
        material
    end

    methods (Access = public)

        function obj = LHSIntegratorStokes(cParams)
            obj.init(cParams);
        end

        function LHS = compute(obj)
            velLHS = obj.computeVelocityLHS();
            D      = obj.computeWeakDivergenceMatrix();
            prsLHS = obj.computePressureLHS(D);
            LHS = [velLHS, D; D',prsLHS];
        end

    end

    methods (Access = private)
    
        function init(obj, cParams)
            obj.dt          = cParams.dt;
            obj.mesh        = cParams.mesh;
            obj.material    = cParams.material;
            obj.pressureFun = cParams.pressureFun;
            obj.velocityFun = cParams.velocityFun;
        end

        function LHS = computeVelocityLHS(obj)
            K = obj.computeVelocityLaplacian();
            M = obj.computeMassMatrix();
            LHS = K + M;
        end

        function D = computeWeakDivergenceMatrix(obj)
            f = @(p,v) DP(v,Grad(p));
            D = IntegrateLHS(f,obj.velocityFun,obj.pressureFun,obj.mesh);            
        end

        function BB = computePressureLHS(obj,D)
            sz = size(D, 2);
            BB = sparse(sz,sz);
        end

        function lhs = computeVelocityLaplacian(obj)
            f = @(u,v) DDP(Grad(v),Grad(u));
            lhs = IntegrateLHS(f,obj.velocityFun,obj.velocityFun,obj.mesh);
            lhs = 1/2*(lhs+lhs');
        end

        function M = computeMassMatrix(obj)
            f = @(u,v) DP(v,u)./obj.dt;
            M = IntegrateLHS(f,obj.velocityFun,obj.velocityFun,obj.mesh,3);
        end

    end

end