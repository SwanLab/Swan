classdef ThreeFieldsComputer < handle

    properties (Access = private)
        mesh
        singleBdMesh
        uFun
        lambdaFun
        material
        dirichletFun
        uGamma
    end

    methods (Access = public)

        function obj = ThreeFieldsComputer(cParams)
            obj.init(cParams)
        end

        function [u,lambda,K,Kc] = solve(obj)
            K   = obj.computeKfine();
            A   = obj.computeConditionMatrix(obj.lambdaFun,obj.uFun);
            LHS = obj.computeLHS(K);
            RHS = obj.computeRHS();
            sol = LHS\RHS;
            [u,lambda] = obj.computeFunctions(sol);
            Kc = u.'*K*u;
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.mesh         = cParams.mesh;
            obj.singleBdMesh = obj.mesh.createSingleBoundaryMesh();
            obj.uFun         = cParams.uFun;
            obj.lambdaFun    = cParams.lambdaFun;
            obj.material     = cParams.material;
            % obj.dirichletFun = cParams.dirichletFun;
            obj.uGamma       = cParams.uGamma;
        end

        function K = computeKfine(obj)
            C = obj.material;
            f = @(u,v) DDP(SymGrad(v),DDP(C,SymGrad(u)));
            K = IntegrateLHS(f,obj.uFun,obj.uFun,obj.mesh,'Domain',2);
        end

        function A = computeLambdaUdom(obj)
             % bMesh = obj.lambdaFun.mesh;
             % uBdFun = obj.uFun.restrictToBoundary();
             f = @(v,u) DP(v,u);
             A= IntegrateLHS(f,obj.lambdaFun,obj.uFun,obj.singleBdMesh,'Boundary',2);
        end

        function A = computeConditionMatrix(obj,u,v)
             f = @(v,u) DP(v,u);
             A= IntegrateLHS(f,v,u,obj.singleBdMesh,'Boundary',2);
        end

        function LHS=computeLHS(obj,K)
            f = @(u,v) DP(v,u);
            C = IntegrateLHS(f,obj.uFun,obj.lambdaFun,obj.mesh,'Boundary',2); 
            Z  = zeros(obj.lambdaFun.nDofs);
            LHS = [K C; C' Z];
        end   

        function RHS = computeRHS(obj)            
            uD   = obj.dirichletFun;
            rDir = zeros(obj.lambdaFun.nDofs,numel(uD));
            for iD = 1:numel(uD)
                f = @(v) DP(v,uD{iD});
                rDir(:,iD)= IntegrateRHS(f,obj.lambdaFun,obj.mesh,'Boundary',2);
            end
            Z   = zeros(obj.uFun.nDofs,numel(uD));
            RHS = [Z; rDir];
        end      

        function [uFun,lambdaFun] = computeFunctions(obj,sol)
            u = sol(1:obj.uFun.nDofs,:);
            L = -sol(obj.uFun.nDofs+1:end,:);
            u=full(u);
            L=full(L);
            %loop
            % u = copy(obj.uFun);
            % uFun{i} = u.setFValues();
            %u, L as uFun and lambdFun

        end

    end

end