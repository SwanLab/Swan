classdef ElasticHarmonicExtension < handle

    properties (Access = private)
        mesh
        boundaryMesh
        uFun
        lambdaFun
        material
        dirichletFun
        localGlobalConnecBd
    end

    methods (Access = public)

        function obj = ElasticHarmonicExtension(cParams)
            obj.init(cParams)

        end

        function [u,L,K] = solve(obj)
            K = obj.computeKfine();
            LHS=obj.computeLHS(K);
            RHS=obj.computeRHS();
            sol = LHS\RHS;
            u = sol(1:obj.uFun.nDofs,:);
            L = -sol(obj.uFun.nDofs+1:end,:);
            if isa(obj.lambdaFun, "LagrangianFunction")
                l2g_dof = ((obj.localGlobalConnecBd*obj.uFun.ndimf)' - ((obj.uFun.ndimf-1):-1:0))';
                l2g_dof = l2g_dof(:);
                uB = u(l2g_dof, :);
                L = uB'*L;
            end
            u=full(u);
            L=full(L);
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.mesh         = cParams.mesh;
            obj.boundaryMesh = cParams.boundaryMesh;
            obj.uFun         = cParams.uFun;
            obj.lambdaFun    = cParams.lambdaFun;
            obj.material     = cParams.material;
            obj.dirichletFun = cParams.dirichletFun;
            obj.localGlobalConnecBd = cParams.localGlobalConnecBd;
        end

        function K = computeKfine(obj)
            C = obj.material;
            f = @(u,v) DDP(SymGrad(v),DDP(C,SymGrad(u)));
            K = IntegrateLHS(f,obj.uFun,obj.uFun,obj.mesh,'Domain',2);
        end

        function LHS=computeLHS(obj,K)
            f = @(u,v) DP(v,u);
            C = IntegrateLHS(f,obj.uFun,obj.lambdaFun,obj.mesh,'Boundary',2); 
            Z  = zeros(obj.lambdaFun.nDofs);
            LHS = [K C; C' Z];
        end   

        function RHS=computeRHS(obj)
            test   = LagrangianFunction.create(obj.boundaryMesh, obj.mesh.ndim, 'P1');
            rDir = [];
            uD = obj.dirichletFun;
            for i=1:numel(uD)
                f = @(v) DP(v,uD{i});
                rDire = IntegrateRHS(f,test,obj.boundaryMesh,'Domain',2);
                rDir = [rDir rDire];
            end
            Z   = zeros(obj.uFun.nDofs, size(rDir, 2));
            RHS = [Z; rDir];
        end        

    end

end