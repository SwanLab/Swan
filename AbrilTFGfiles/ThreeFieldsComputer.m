classdef ThreeFieldsComputer < handle

    properties (Access = private)
        mesh
        singleBdMesh
        uFun
        lambdaFun
        material
        %dirichletFun
        uGamma
    end

    methods (Access = public)

        function obj = ThreeFieldsComputer(cParams)
            obj.init(cParams)
        end

        function [u,lambda,K,Kc] = solve(obj)
            K   = obj.computeKfine(); % Matriu de rigidesa (domini)
            A   = obj.computeConditionMatrix(obj.lambdaFun,obj.uFun);% Acomplament de u-lambda (frontera)
            L   = obj.computeConditionMatrix(obj.uGamma,obj.lambdaFun); % Acoplament de lambda-uGamma (Frontera)
            LHS = obj.computeLHS(K,A);
            RHS = obj.computeRHS(L);
            sol = LHS\RHS;
            [u,lambda] = obj.computeFunctions(sol);
            Kc = u.'*K*u; % Que es la Kcoarse?
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.mesh         = cParams.mesh;
            obj.singleBdMesh = obj.mesh.createSingleBoundaryMesh();
            obj.uFun         = cParams.uFun;
            obj.lambdaFun    = cParams.lambdaFun;
            obj.material     = cParams.material;
            obj.uGamma       = cParams.uGamma;
          % obj.dirichletFun = cParams.dirichletFun;

        end

        function K = computeKfine(obj)
            C = obj.material;
            f = @(u,v) DDP(SymGrad(v),DDP(C,SymGrad(u)));
            K = IntegrateLHS(f,obj.uFun,obj.uFun,obj.mesh,'Domain',2);
        end

        function M = computeConditionMatrix(obj,u,v)
             f = @(v,u) DP(v,u);
             M = IntegrateLHS(f,v,u,obj.singleBdMesh,'Boundary',2);
        end

        function LHS = computeLHS(obj,K,A)
            nLambda = size(A, 2);
            Z       = sparse(nLambda, nLambda);
            LHS     = [K A; A' Z];
            %Z  = zeros(obj.lambdaFun.nDofs);
            %f = @(u,v) DP(v,u);
            %C = IntegrateLHS(f,obj.uFun,obj.lambdaFun,obj.mesh,'Boundary',2); 
        end

        function RHS = computeRHS(obj, L)
            nU      = obj.uFun.nDofs;
            nUGamma = size(L,2);
            Z       = sparse(nU, nUGamma);
            RHS     = [Z;L];
            % uD   = obj.dirichletFun;
            % rDir = zeros(obj.lambdaFun.nDofs,numel(uD));
            % for iD = 1:numel(uD)
            %     f = @(v) DP(v,uD{iD});
            %     rDir(:,iD)= IntegrateRHS(f,obj.lambdaFun,obj.mesh,'Boundary',2);
            % end
            % Z   = zeros(obj.uFun.nDofs,numel(uD));
            % RHS = [Z; L];
        end  
                            
        function [u,lambda] = computeFunctions(obj,sol)
            nU     = obj.uFun.nDofs;
            u      = sol(1:nU, :);
            lambda = -sol(nU+1:end, :);
            u      = full(u);
            lambda = full(lambda); 
            %loop
            % u = copy(obj.uFun);
            % uFun{i} = u.setFValues();
            %u, L as uFun and lambdFun
        end

        % function A = computeLambdaUdom(obj)
        %      % bMesh = obj.lambdaFun.mesh;
        %      % uBdFun = obj.uFun.restrictToBoundary();
        %      f = @(v,u) DP(v,u);
        %      A= IntegrateLHS(f,obj.lambdaFun,obj.uFun,obj.singleBdMesh,'Boundary',2);
        % end

        % function RHS = computeRHS(obj)            
        %     uD   = obj.dirichletFun;
        %     rDir = zeros(obj.lambdaFun.nDofs,numel(uD));
        %     for iD = 1:numel(uD)
        %         f = @(v) DP(v,uD{iD});
        %         rDir(:,iD)= IntegrateRHS(f,obj.lambdaFun,obj.mesh,'Boundary',2);
        %     end
        %     Z   = zeros(obj.uFun.nDofs,numel(uD));
        %     RHS = [Z; rDir];
        % end  

        % function [uFun,lambdaFun] = computeFunctions(obj,sol)
        %     u = sol(1:obj.uFun.nDofs,:);
        %     L = -sol(obj.uFun.nDofs+1:end,:);
        %     u=full(u);
        %     L=full(L);
        %     %loop
        %     % u = copy(obj.uFun);
        %     % uFun{i} = u.setFValues();
        %     %u, L as uFun and lambdFun
        % end

    end

end