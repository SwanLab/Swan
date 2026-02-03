classdef ThreeFieldsComputer2 < handle
    properties (Access = private)
        mesh1, mesh2
        bdSubmesh
        uFun, uFun1, uFun2
        uFun1bd, uFun2bd
        lambdaFun1, lambdaFun2
        material1, material2
        uGamma1, uGamma2
        bc1, bc2
        localGlobalConnec       
        fixedDOFs 
    end
    methods (Access = public)
        function obj = ThreeFieldsComputer2(cParams)
            obj.init(cParams)
        end
        function [u, lambda1, lambda2, K, Kc] = solve(obj)
            K1 = obj.computeKfine(obj.uFun1, obj.material1);
            K2 = obj.computeKfine(obj.uFun2, obj.material2);
            
            A1 = obj.computeConditionMatrix(obj.lambdaFun1, obj.uFun1bd, obj.uFun1, obj.bdSubmesh{1,1}{2});
            A2 = obj.computeConditionMatrix(obj.lambdaFun2, obj.uFun2bd, obj.uFun2, obj.bdSubmesh{1,2}{1});
                       
            L1 = obj.computeConditionMatrix2(obj.lambdaFun1, obj.uGamma1);
            L2 = obj.computeConditionMatrix2(obj.lambdaFun2, obj.uGamma2);
            
            LHS = obj.computeLHS(K1, K2, A1, A2, L1, L2);
            RHS = obj.computeRHS(L1, L2);

            RHS = obj.applyNeumann(RHS);

            [LHS, RHS] = obj.applyDirichlet(LHS, RHS);            
                       
            sol = LHS \ RHS;
            [u, lambda1, lambda2] = obj.computeFunctions(sol);
            
            K = blkdiag(K1, K2);
            Kc = u.' * K * u;
        end
    end
    methods (Access = private)
        function init(obj, cParams)
            obj.mesh1 = cParams.mesh1;
            obj.mesh2 = cParams.mesh2;
            obj.bdSubmesh = cParams.bdSubmesh;
            obj.uFun = cParams.uFun;
            obj.uFun1 = cParams.uFun1;
            obj.uFun2 = cParams.uFun2;
            obj.uFun1bd = cParams.uFun1bd;
            obj.uFun2bd = cParams.uFun2bd;
            obj.lambdaFun1 = cParams.lambdaFun1;
            obj.lambdaFun2 = cParams.lambdaFun2;
            obj.material1 = cParams.material1;
            obj.material2 = cParams.material2;
            obj.uGamma1 = cParams.uGamma1;
            obj.uGamma2 = cParams.uGamma2; 
            obj.bc1 = cParams.bc1;
            obj.bc2 = cParams.bc2;
        end
        
        function K = computeKfine(obj, uFun, material)
            C = material;
            f = @(u,v) DDP(SymGrad(v), DDP(C, SymGrad(u)));
            K = IntegrateLHS(f, uFun, uFun, uFun.mesh, 'Domain', 2);
        end

        function Mglob = computeConditionMatrix(obj, testFun, trialFun, uDomain, bdSubmesh) 
            f = @(v,u) DP(v,u);
            M = IntegrateLHS(f, testFun, trialFun, testFun.mesh, 'Domain', 2);
        
            nDim        = testFun.ndimf;
            localNodes  = unique(testFun.mesh.connec(:), 'Stable');
            globalNodes = bdSubmesh.globalConnec(:);
        
            % Expand nodes to DOFs: node k -> DOFs [(k-1)*nDim+1, ..., k*nDim]
            localDofs  = (localNodes - 1) * nDim + (1:nDim);   % each row: DOFs of one node
            globalDofs = (globalNodes - 1) * nDim + (1:nDim);
        
            localDofs  = localDofs(:);   % flatten, preserving node-major order
            globalDofs = globalDofs(:);
        
            % Assemble: rows stay local (Lambda), cols map to global (U)
            Mglob = sparse(testFun.nDofs, uDomain.nDofs);
            Mglob(:, globalDofs) = M(:, localDofs);
        end
        
        function M2 = computeConditionMatrix2(obj, testFun, trialFun)
            f = @(v,u) DP(v,u);
            integrationMesh = testFun.mesh; 
            M2 = IntegrateLHS(f, testFun, trialFun, integrationMesh, 'Domain', 2);

        end
        
        function LHS = computeLHS(obj, K1, K2, A1, A2, L1, L2)
            nU1  = size(K1, 1);
            nU2  = size(K2, 1);
            nL1  = size(A1, 1);
            nL2  = size(A2, 1);
            nUG  = size(L1, 2);   % L1 is (nL1 x nUG)
        
            LHS = [K1,               sparse(nU1,nU2), A1',              sparse(nU1,nL2), sparse(nU1,nUG);
                   sparse(nU2,nU1),  K2,              sparse(nU2,nL1),  A2',            sparse(nU2,nUG);
                   A1,               sparse(nL1,nU2), sparse(nL1,nL1),  sparse(nL1,nL2),-L1;
                   sparse(nL2,nU1),  A2,              sparse(nL2,nL1),  sparse(nL2,nL2),-L2;
                   sparse(nUG,nU1),  sparse(nUG,nU2), -L1',             -L2',           sparse(nUG,nUG)];
        end
        
        function RHS = computeRHS(obj, L1, L2)
            nU1  = obj.uFun1.nDofs;
            nU2  = obj.uFun2.nDofs;
            nL1  = obj.lambdaFun1.nDofs;
            nL2  = obj.lambdaFun2.nDofs;
            nUG  = obj.uGamma1.nDofs;
            totalDofs = nU1 + nU2 + nL1 + nL2 + nUG;
            nModes    = 8;
            RHS       = sparse(totalDofs, nModes);
            startIdx  = nU1 + nU2 + nL1 + nL2 + 1;
            for i = 1:min(nModes, nUG)
                RHS(startIdx + i - 1, i) = 1.0;
            end
        end

        function RHS = applyNeumann(obj, RHS)
            % Neumann force lives in the subdomain 2 block
            % Offset = nDofs of subdomain 1
            nU1   = obj.uFun1.nDofs;
            t     = obj.bc2.tractionFun;
            rhs2  = zeros(obj.uFun2.nDofs, 1);
            if ~isempty(t)
                for i = 1:numel(t)
                    rhs2 = rhs2 + t(i).computeRHS(obj.uFun2);
                end
            end
            % Add to all nModes columns
            RHS(nU1+1:nU1+obj.uFun2.nDofs, :) = RHS(nU1+1:nU1+obj.uFun2.nDofs, :) + rhs2;
        end
        
        function [LHS, RHS] = applyDirichlet(obj, LHS, RHS)
            dirichDofs = obj.bc1.dirichlet_dofs;
            dirichVals = obj.bc1.dirichlet_vals;
   
            if isempty(dirichDofs)
                return;
            end
        
            nModes = size(RHS, 2);
            RHS = RHS - LHS(:, dirichDofs) * dirichVals * ones(1, nModes);
        
            LHS(dirichDofs, :) = 0;
            LHS(:, dirichDofs) = 0;
            LHS(sub2ind(size(LHS), dirichDofs, dirichDofs)) = 1;
        
            RHS(dirichDofs, :) = dirichVals * ones(1, nModes);
        end
        

        function [u, lambda1, lambda2] = computeFunctions(obj, sol)
            nU1 = obj.uFun1.nDofs;
            nU2 = obj.uFun2.nDofs;
            nL1 = obj.lambdaFun1.nDofs;
            u1      = sol(1:nU1, :);
            u2      = sol(nU1+1:nU1+nU2, :);
            lambda1 = sol(nU1+nU2+1:nU1+nU2+nL1, :);
            lambda2 = sol(nU1+nU2+nL1+1:nU1+nU2+nL1+obj.lambdaFun2.nDofs, :);
            u       = full([u1; u2]);
            lambda1 = full(lambda1);
            lambda2 = full(lambda2);
        end
    end
end