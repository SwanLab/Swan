classdef ElasticHarmonicExtension < handle

    properties (Access = private)
        mesh
        uFun
        lambdaFun
        material
        dirichletFun
        bMesh 
        type
        nLambda
        lambDofsGlobal
    end

    methods (Access = public)

        function obj = ElasticHarmonicExtension(cParams)
            obj.init(cParams);
        end

        function [u, lambda, K, Kc] = solve(obj)
            RHS = obj.computeRHS();                       
            [LHS,K] = obj.computeLHS();          
            [u, lambda] = obj.computeSolution(LHS,RHS);
            Kc  = u.' * K * u;
        end

    end

    methods (Access = private)

        function init(obj, cParams)
            obj.mesh         = cParams.mesh;
            obj.uFun         = cParams.uFun;
            obj.material     = cParams.material;
            obj.dirichletFun = cParams.dirichletFun;
            obj.type         = cParams.type;

            % Dynamic strategy dispatch for type initialization
            methodName = sprintf('init%s', obj.type);
            obj.(methodName)(cParams);
        end

        function initContinuous(obj, cParams)
            obj.lambdaFun = cParams.lambdaFun;
        end

        function initDiscontinuous(obj, cParams)
            if ~iscell(cParams.lambdaFun)
                obj.lambdaFun     = {cParams.lambdaFun};
                obj.bMesh{1}.mesh = cParams.bMesh;
            else
                obj.lambdaFun = cParams.lambdaFun;
                obj.bMesh     = cParams.bMesh;
            end
            offset = 0;
            for i = 1:numel(obj.lambdaFun)
                nLi = obj.lambdaFun{i}.nDofs;
                obj.lambDofsGlobal{i} = offset + (1:nLi);
                offset = offset + nLi;
            end
            obj.nLambda = offset;
        end

        function K = computeKfine(obj)
            C = obj.material;
            f = @(u, v) DDP(SymGrad(v), DDP(C, SymGrad(u)));
            K = IntegrateLHS(f, obj.uFun, obj.uFun, obj.mesh, 'Domain', 2);
        end

        function [LHS,K] = computeLHS(obj)
            K          = obj.computeKfine();
            methodName = sprintf('computeLHS%s', obj.type);
            LHS        = obj.(methodName)(K);
        end   

        function LHS = computeLHSContinuous(obj, K)
            f = @(u, v) DP(v, u);
            C = IntegrateLHS(f, obj.uFun, obj.lambdaFun, obj.mesh, 'Boundary', 2);
            Z = sparse(obj.lambdaFun.nDofs, obj.lambdaFun.nDofs);
            LHS = [K C; C' Z];
        end

        function LHS = computeLHSDiscontinuous(obj, K)
            f = @(u, v) DP(v, u);
            C = sparse(obj.uFun.nDofs, obj.nLambda);
            [uBfun, iGlob] = obj.restrictToBoundary(obj.uFun);

            for i = 1:numel(obj.lambdaFun)
                jGlob = obj.lambDofsGlobal{i};
                lhsLoc = IntegrateLHS(f, uBfun{i}, obj.lambdaFun{i}, obj.bMesh{i}.mesh, 'Domain', 2);
                [iLoc, jLoc, vals] = find(lhsLoc);
                Ci = sparse(iGlob{i}(iLoc), jGlob(jLoc), vals, obj.uFun.nDofs, obj.nLambda);
                C  = C + Ci;
            end
            Z   = sparse(obj.nLambda, obj.nLambda);
            LHS = [K C; C' Z];
        end

        function RHS = computeRHS(obj)   
            methodName = sprintf('computeRHS%s', obj.type);
            rDir = obj.(methodName)();
            Z    = sparse(obj.uFun.nDofs, numel(obj.dirichletFun));
            RHS  = [Z; rDir];
        end      

        function rDir = computeRHSContinuous(obj)
            uD   = obj.dirichletFun;
            rDir = zeros(obj.lambdaFun.nDofs, numel(uD));
            for iD = 1:numel(uD)
                f = @(v) DP(v, uD{iD});
                rDir(:, iD) = IntegrateRHS(f, obj.lambdaFun, obj.mesh, 'Boundary', 2);
            end
        end

        function rDir = computeRHSDiscontinuous(obj)
            uD = obj.dirichletFun;
            nBoundaryEntities  = numel(obj.lambdaFun);
            nCoarsePerBoundary = numel(uD) / nBoundaryEntities;

            iLambda0 = 0;
            rDir = zeros(obj.nLambda, numel(uD));

            for iBoundary = 1:nBoundaryEntities
                nLi = obj.lambdaFun{iBoundary}.nDofs;
                coarseIds = (iBoundary - 1) * nCoarsePerBoundary + (1:nCoarsePerBoundary);

                for iD = coarseIds
                    f  = @(v) DP(v, uD{iD});
                    ri = IntegrateRHS(f, obj.lambdaFun{iBoundary}, obj.bMesh{iBoundary}.mesh, 'Domain', 2);
                    rDir(iLambda0 + (1:nLi), iD) = ri;
                end
                iLambda0 = iLambda0 + nLi;
            end
        end

        function [u, lambda] = computeSolution(obj, LHS,RHS)
            sol = LHS \ RHS;
            u      = sol(1:obj.uFun.nDofs, :);
            lambda = -sol(obj.uFun.nDofs+1:end, :);
        end

        function [bFun, gFunc] = restrictToBoundary(obj, fun)
            bFun  = cell(1, numel(obj.bMesh));
            gFunc = bFun;
            for i = 1:numel(obj.bMesh)
                l2g      = unique(obj.bMesh{i}.globalConnec(:));
                lastDofs = (l2g * fun.ndimf)';
                l2gdof   = zeros(length(lastDofs), fun.ndimf);
                for j = 1:fun.ndimf
                    l2gdof(:, j) = lastDofs - (fun.ndimf - j);
                end
                bFun{i} = LagrangianFunction.create(obj.bMesh{i}.mesh, fun.ndimf, fun.order);
                val = fun.fValues(l2gdof);
                bFun{i}.setFValues(val);
                l2g_map  = reshape(l2gdof', [], 1);
                gFunc{i} = @(iLoc) l2g_map(iLoc);
            end
        end

    end

end