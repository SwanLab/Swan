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
    end

    methods (Access = public)

        function obj = ElasticHarmonicExtension(cParams)
            obj.init(cParams)
        end

        function [u,lambda,K,Kc] = solve(obj)
            RHS = obj.computeRHS();            
            K   = obj.computeKfine();
            LHS = obj.computeLHS(K);
            sol = LHS\RHS;
            [u,lambda] = obj.computeFunctions(sol);
            Kc = u.'*K*u;
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.mesh           = cParams.mesh;
            % obj.lambdaFun    = cParams.lambdaFun;
            obj.uFun           = cParams.uFun;
            obj.material       = cParams.material;
            obj.dirichletFun   = cParams.dirichletFun;
            obj.type           = cParams.type;

            switch obj.type
                case 'continuous'
                    obj.lambdaFun  = cParams.lambdaFun;
                    % obj.bMesh      = cParams.bMesh;
                case 'discontinuous'
                    if ~iscell(cParams.lambdaFun)
                        obj.lambdaFun      = {cParams.lambdaFun};
                        obj.bMesh{1}.mesh  = cParams.bMesh;
                    else
                        obj.lambdaFun  = cParams.lambdaFun;
                        obj.bMesh      = cParams.bMesh;
                    end
                   obj.nLambda = 0;
                    for iLam = 1:numel(obj.lambdaFun)
                        obj.nLambda = obj.nLambda + obj.lambdaFun{iLam}.nDofs;
                    end
            end            
        end


        function K = computeKfine(obj)
            C = obj.material;
            f = @(u,v) DDP(SymGrad(v),DDP(C,SymGrad(u)));
            K = IntegrateLHS(f,obj.uFun,obj.uFun,obj.mesh,'Domain',2);
        end

        function LHS=computeLHS(obj,K)
            f = @(u,v) DP(v,u);

            switch obj.type
                case 'continuous'
                    C = IntegrateLHS(f,obj.uFun,obj.lambdaFun,obj.mesh,'Boundary',2);
                    Z  = sparse(obj.lambdaFun.nDofs,obj.lambdaFun.nDofs);
                case 'discontinuous'
                    C= sparse(obj.uFun.nDofs,obj.nLambda);
                    [uBfun,iGlob]=obj.restrictToBoundary(obj.uFun);
                    jGlob= @(iLoc) iLoc;
                    for i = 1:numel(obj.lambdaFun)
                        lhsLoc = IntegrateLHS(f,uBfun{i},obj.lambdaFun{i},obj.bMesh{i}.mesh,'Domain',2);
                        [iLoc,jLoc,vals] = find(lhsLoc);
                        Ci = sparse(iGlob{i}(iLoc),jGlob(jLoc),vals, obj.uFun.nDofs,obj.nLambda);
                        C=C+Ci;
                    end
                    Z  = sparse(obj.nLambda,obj.nLambda);
            end
            LHS = [K C; C' Z];
        end   

        function RHS = computeRHS(obj)   
            uD   = obj.dirichletFun;
            
            switch obj.type
                case 'continuous'
                    rDir = zeros(obj.lambdaFun.nDofs,numel(uD));
                    for iD = 1:numel(uD)
                        f = @(v) DP(v,uD{iD});
                        rDir(:,iD)= IntegrateRHS(f,obj.lambdaFun,obj.mesh,'Boundary',2);
                    end

                case 'discontinuous'
                    nCoarsePerSegment = 4;
                    iLambda0 = 0;
                    rDir = zeros(obj.nLambda,numel(uD));

                    for iSegment = 1:numel(obj.lambdaFun)
                        nLi = obj.lambdaFun{iSegment}.nDofs;
                        coarseIds = (iSegment-1)*nCoarsePerSegment + (1:nCoarsePerSegment);

                        for iD = coarseIds
                            f = @(v) DP(v,uD{iD});
                            ri = IntegrateRHS(f,obj.lambdaFun{iSegment},obj.bMesh{iSegment}.mesh,'Domain',2);
                            rDir(iLambda0 + (1:nLi),iD) = ri;
                        end
                        iLambda0 = iLambda0 + nLi;
                    end
                    
                    % for iD = 1:numel(uD)
                    %     n = 0;
                    %     iSegment = floor((iD-1)/4)+1;
                    % 
                    %     for iLam = 1:numel(obj.lambdaFun)
                    %         nLi = obj.lambdaFun{iLam}.nDofs;
                    %         if iLam == iSegment
                    %             f = @(v) DP(v,uD{iD});
                    %             ri = IntegrateRHS(f,obj.lambdaFun{iLam},obj.bMesh{iLam}.mesh,'Domain',2);
                    %             rDir(n+(1:nLi),iD) = ri;
                    %         end
                    %         n = n + nLi;
                    %     end
                    % end
            end

            Z   = sparse(obj.uFun.nDofs,numel(uD));
            RHS = [Z; rDir];
        end      

        function [uFun,lambdaFun] = computeFunctions(obj,sol)
            u = sol(1:obj.uFun.nDofs,:);
            L = -sol(obj.uFun.nDofs+1:end,:);
            u=full(u);
            L=full(L);

            % uFun=cell(1,size(u,2));
            % lambdaFun=cell(1,size(L,2));
            % 
            % for i=1:size(u,2)
            %     uFun{i}=copy(obj.uFun);
            %     fV= reshape(u(:,i),[obj.mesh.ndim,obj.mesh.nnodes])';
            %     uFun{i}.setFValues(fV);
            % end
            % 
            % for i=1:size(L,2)
            %     lambdaFun{i}=copy(obj.lambdaFun);
            %     fV= reshape(L(:,i),obj.mesh.ndim,[]).';
            %     lambdaFun{i}.setFValues(fV);
            % end
            
            uFun=u;
            lambdaFun=L;

        end

        function uFun = createDispFun(obj)
            ntest  = size(obj.dirichletFun,2);
            s.mesh  = obj.mesh;
            s.ndimf = obj.mesh.ndim;
            s.order = 'P1';
            for i = 1:ntest
                fValues   = obj.fValuesTraining(:,i);
                s.fValues = reshape(fValues,obj.mesh.ndim,[])' ;
                uFun(i)   = LagrangianFunction(s);
            end
        end

        function [bFun,gFunc]=restrictToBoundary(obj,fun)
            bFun=cell(1,numel(obj.bMesh));
            gFunc=bFun;
            for i=1:numel(obj.bMesh)
                l2g= unique(obj.bMesh{i}.globalConnec(:));
                lastDofs = (l2g * fun.ndimf)';
                l2gdof = zeros(length(lastDofs),fun.ndimf);
                for j = 1:fun.ndimf
                    l2gdof(:,j) = lastDofs - (fun.ndimf-j);
                end
                bFun{i} = LagrangianFunction.create(obj.bMesh{i}.mesh,fun.ndimf,fun.order);
                val = fun.fValues(l2gdof);
                bFun{i}.setFValues(val);
                l2g_map = reshape(l2gdof',[],1);
                gFunc{i} = @(iLoc) l2g_map(iLoc);
            end
        end
    end

end