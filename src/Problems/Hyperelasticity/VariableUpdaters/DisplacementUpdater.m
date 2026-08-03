classdef DisplacementUpdater < handle
    
    properties (Access = private)
        functional
        monitor

        tol
        maxIter

        pcgIterHistoryThisStep
        pcgMeanHistoryPerStep

        %EIFEM / ILU Constants
        eifemData
        activePreconditioner
        compareEIFEM
        useConstantPreconditioners
        MeifemConstant
        MiluConstant
        LHSforILUConstant

    end

    methods (Access = public)

        function obj = DisplacementUpdater(cParams)
            obj.init(cParams);
        end

        function [u,rFun,costArray,iter] = update(obj,u,bc,costArray)

            obj.pcgIterHistoryThisStep = [];

            i = 0; err = 1; costOld = costArray(end);
            while (abs(err) > obj.tol) && (i < obj.maxIter)
                LHS = obj.functional.computeHessian(u);
                RHS = obj.functional.computeGradient(u,bc);
                u.setFValues(obj.computeDisplacement(LHS,RHS,u,bc));

                [err, cost] = obj.computeErrorCost(u,bc,costOld);
                costArray(end+1) = cost;
                costOld = cost;

                i = i+1;
                obj.monitor.printCost('iterU',i,cost,err);
                obj.monitor.update(length(costArray),{[],[cost],[],[]});
                % obj.monitor.refresh(); 
            end
         
            if ~isempty(obj.pcgIterHistoryThisStep)
                meanCG_thisStep = mean(obj.pcgIterHistoryThisStep);
            else               
                meanCG_thisStep = NaN;
            end

            obj.pcgMeanHistoryPerStep(end+1) = meanCG_thisStep;

            fprintf('Load step acabat: Newton iters = %d, PCG mean iters/Newton = %.2f\n', i, meanCG_thisStep);

            figure(202); clf;
            bar(obj.pcgMeanHistoryPerStep); hold on; grid on;
            xlabel('Load step');
            ylabel('Mean PCG iterations per Newton');
            
            overallMean = mean(obj.pcgMeanHistoryPerStep,'omitnan');
            yline(overallMean, 'r', 'LineWidth', 2);

            title(sprintf('Mean PCG iters per Newton, per load step (overall mean = %.2f)', overallMean));
            legend('Mean per load step','Overall mean','Location','best');  

            rFun = obj.computeReactions(LHS,u,bc);
            iter = i;
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.functional = cParams.functional;
            obj.monitor    = cParams.monitor;
            obj.tol        = cParams.tolerance;
            obj.maxIter    = cParams.maxIter;
            obj.pcgIterHistoryThisStep = [];
            obj.pcgMeanHistoryPerStep  = [];

            % %EIFEM
            % obj.eifemData = cParams.eifemData;
            % obj.activePreconditioner = cParams.activePreconditioner;
            % obj.compareEIFEM = cParams.compareEIFEM;

            % EIFEM / ILU constants
            obj.eifemData = cParams.eifemData;
            obj.activePreconditioner = cParams.activePreconditioner;
            obj.compareEIFEM = cParams.compareEIFEM;
            
            if isfield(cParams,'useConstantPreconditioners')
                obj.useConstantPreconditioners = cParams.useConstantPreconditioners;
            else
                obj.useConstantPreconditioners = false;
            end

            % Inicialització precondicionadors constants
            obj.MeifemConstant = [];
            obj.MiluConstant   = [];
            obj.LHSforILUConstant = [];

        end

        function uOut = computeDisplacement(obj,LHSfull, RHSfull,uIn,bc)
            [LHS,RHS] = fullToReduced(obj,LHSfull,RHSfull,bc);
            if ~isempty(LHS)
                uInVec = reshape(uIn.fValues',[uIn.nDofs 1]);
                uOutVec = uInVec;

                uInFree = uInVec(bc.free_dofs);
                uOutFree = obj.updateWithNewton(LHS,RHS,uInFree,bc);
                uOutVec(bc.free_dofs) = uOutFree;
                uOut = reshape(uOutVec,[flip(size(uIn.fValues))])';
            else
                uOut = uIn.fValues;
            end
        end

        function [LHS,RHS] = fullToReduced(~,LHS,RHS,bc)
            free_dofs = bc.free_dofs;
            LHS = LHS(free_dofs, free_dofs);
            RHS = RHS(free_dofs);
        end

        function xNew = updateWithNewton(obj,LHS,RHS,x,bc)

            % disp(obj.compareEIFEM)                  % debug
            % disp(obj.useConstantPreconditioners)    % debug

            tol = 1e-4;
            maxIter = 100;

            b  = -RHS;
            x0 = zeros(size(b));

            % Solver directe (referència)
            t = tic;
            deltaX_direct = LHS \ b;
            timeDirect = toc(t);

            % Precondicionador ILU variable
            % De moment el desactivem perquè volem precondicionadors constants
            % Milu = obj.createILUpreconditioner(LHS);

            results = struct([]);

            % 1) PCG
            t = tic;
            [dx,flag,relres,iter,resvec] = pcg(LHS,b,tol,maxIter,[],[],x0);
            results(1).name   = 'PCG';
            results(1).dx     = dx;
            results(1).flag   = flag;
            results(1).relres = relres;
            results(1).iter   = iter;
            results(1).resvec = resvec;
            results(1).time   = toc(t);

            % 2) PCG + ILU
            % t = tic;
            % [dx,flag,relres,iter,resvec] = pcg(LHS,b,tol,maxIter,Milu,[],x0);
            % results(2).name   = 'PCG_ILU';
            % results(2).dx     = dx;
            % results(2).flag   = flag;
            % results(2).relres = relres;
            % results(2).iter   = iter;
            % results(2).resvec = resvec;
            % results(2).time   = toc(t);

            % 2) PCG + ILU constant
            Milu = obj.getConstantILUpreconditioner(LHS);
            
            t = tic;
            [dx,flag,relres,iter,resvec] = pcg(LHS,b,tol,maxIter,Milu,[],x0);
            results(2).name   = 'PCG_ILU_CONSTANT';
            results(2).dx     = dx;
            results(2).flag   = flag;
            results(2).relres = relres;
            results(2).iter   = iter;
            results(2).resvec = resvec;
            results(2).time   = toc(t);
            
            
            if obj.compareEIFEM
                
                if obj.useConstantPreconditioners
                    % EIFEM constant
                    Meifem = obj.getConstantEIFEMpreconditioner(bc);
                else
                    % EIFEM recomputat (com abans)
                    bcApplierCurrent = obj.createCurrentBCApplier(bc);
                    Meifem = obj.createEIFEMpreconditioner(bcApplierCurrent);
                end
            
                t = tic;
                [dx,flag,relres,iter,resvec] = pcg(LHS,b,tol,maxIter,Meifem,[],x0);


                results(end+1).name   = 'PCG_EIFEM';
                results(end).dx       = dx;
                results(end).flag     = flag;
                results(end).relres   = relres;
                results(end).iter     = iter;
                results(end).resvec   = resvec;
                results(end).time     = toc(t);

                % 4) PCG + ILU-EIFEM-ILU
                % LHSf   = @(v) LHS*v;
                % Mcombo = @(r) Preconditioner.multiplePrec(r,LHSf,Milu,Meifem,Milu);
                %
                % t = tic;
                % [dx,flag,relres,iter,resvec] = pcg(LHS,b,tol,maxIter,Mcombo,[],x0);
                % results(4).name   = 'PCG_ILU_EIFEM_ILU';
                % results(4).dx     = dx;
                % results(4).flag   = flag;
                % results(4).relres = relres;
                % results(4).iter   = iter;
                % results(4).resvec = resvec;
                % results(4).time   = toc(t);

                % 4) PCG + ILU-EIFEM-ILU constant
                if obj.useConstantPreconditioners
                
                    
                    Milu   = obj.getConstantILUpreconditioner(LHS); 
                    Meifem = obj.getConstantEIFEMpreconditioner(bc);
                
                    LHSfun = @(v) LHS*v;
                    Mcombo = @(r) Preconditioner.multiplePrec(r,LHSfun,Milu,Meifem,Milu);
                
                    t = tic;
                    [dx,flag,relres,iter,resvec] = pcg(LHS,b,tol,maxIter,Mcombo,[],x0);
                
                    results(end+1).name   = 'PCG_ILU_EIFEM_ILU_CONSTANT';
                    results(end).dx     = dx;
                    results(end).flag   = flag;
                    results(end).relres = relres;
                    results(end).iter   = iter;
                    results(end).resvec = resvec;
                    results(end).time   = toc(t);
                
                end

            end

            % Guarda les iteracions del solver actiu
            switch obj.activePreconditioner
                case 'PCG'
                    idx = find(strcmp({results.name},'PCG'),1);
                case 'PCG_ILU'
                    idx = find(strcmp({results.name},'PCG_ILU'),1);
                case 'PCG_ILU_CONSTANT'
                    idx = find(strcmp({results.name},'PCG_ILU_CONSTANT'),1);
                case 'PCG_EIFEM'
                    idx = find(strcmp({results.name},'PCG_EIFEM'),1);
                case 'PCG_EIFEM_CONSTANT'
                    idx = find(strcmp({results.name},'PCG_EIFEM'),1);
                case 'PCG_ILU_EIFEM_ILU'
                    idx = find(strcmp({results.name},'PCG_ILU_EIFEM_ILU'),1);
                case 'PCG_ILU_EIFEM_ILU_CONSTANT'
                    idx = find(strcmp({results.name},'PCG_ILU_EIFEM_ILU_CONSTANT'),1);
                otherwise
                    idx = find(strcmp({results.name},'PCG'),1);
            end

            if isempty(idx)
                error('Active preconditioner "%s" has not been computed in results.', obj.activePreconditioner);
            end

            deltaX = results(idx).dx;

            % %DEBUG
            % errDirecte = norm(deltaX - deltaX_direct)/norm(deltaX_direct);
            % fprintf('Error vs solver directe: %.6e\n', errDirecte);
            % %

            obj.pcgIterHistoryThisStep(end+1) = results(idx).iter;

            % Figura 101: residuals + error respecte directe
            figure(101); clf; hold on; grid on;
            legends = {};

            for k = 1:numel(results)
                resPcg = results(k).resvec / norm(b);
                resSol = norm(results(k).dx - deltaX_direct) / norm(deltaX_direct);

                semilogy(resPcg,'LineWidth',1.5);
                semilogy(length(resPcg), resSol, 'o', 'MarkerSize',8, 'LineWidth',1.5);

                legends{end+1} = [results(k).name ' residual'];
                legends{end+1} = [results(k).name ' error'];
            end

            xlabel('iteration');
            ylabel('residual / error');
            legend(legends,'Location','best');
            title('Comparativa de solucionadors lineals');

            % Figura 102: iteracions
            figure(102); clf;
            bar(categorical({results.name}), [results.iter]);
            ylabel('PCG iterations');
            title('Comparativa d''iteracions');

            % Figura 103: temps
            figure(103); clf;
            bar(categorical([{'Direct'}, {results.name}]), [timeDirect, [results.time]]);
            ylabel('CPU time (s)');
            title('Comparativa de temps');

            % Actualització de Newton amb el solver escollit
            xNew = x + deltaX;
        end
          

        function [e, cost] = computeErrorCost(obj,u,bc,costOld)
            cost = obj.functional.computeCost(u,bc); % To include extWork
            % e = (cost - costOld)/cost;

            if abs(cost) < 1e-14
                e = 0;
            else
                e = (cost - costOld)/cost;
            end

        end

        function rFun = computeReactions(~,LHS,u,bc)
            constrainedDofs = bc.dirichlet_dofs;
            uVec = reshape(u.fValues',[u.nDofs 1]);
            KR   = LHS(constrainedDofs,:);
            rVec = zeros(size(uVec));
            rVec(constrainedDofs,1) = -KR*uVec;
            
            s.fValues = reshape(rVec,[flip(size(u.fValues))])';
            s.order   = 'P1';
            s.mesh    = u.mesh;
            rFun = LagrangianFunction(s);
        end

        function Milu = createILUpreconditioner(~,LHS)

            s.LHS = sparse(LHS);
            s.type = 'ILU';
            M = Preconditioner.create(s);
            Milu = @(r) M.apply(r);

        end

        function Milu = getConstantILUpreconditioner(obj,LHS)

            if isempty(obj.MiluConstant)
                obj.LHSforILUConstant = LHS;
                obj.MiluConstant = obj.createILUpreconditioner(obj.LHSforILUConstant);
            end
        
            Milu = obj.MiluConstant;
        
        end

        function Meifem = getConstantEIFEMpreconditioner(obj,bc)

            if isempty(obj.MeifemConstant)
                bcApplierCurrent = obj.createCurrentBCApplier(bc);
                obj.MeifemConstant = obj.createEIFEMpreconditioner(bcApplierCurrent);
            end
        
            Meifem = obj.MeifemConstant;
        
        end


        function Meifem = createEIFEMpreconditioner(obj,bcApplierCurrent)

            d = obj.eifemData;

            s.RVE         = TrainedRVE(d.fileNameEIFEM);
            s.mesh        = obj.createCoarseMesh(d.referenceMesh,d.nSubdomains,d.tolSameNode);
            s.DirCond     = d.dir;
            s.nSubdomains = d.nSubdomains;
            
            % CANVI 4: Reconstruir DirCond per a la mesh coarse
            minx = min(s.mesh.coord(:,1));
            maxx = max(s.mesh.coord(:,1));
            tolSameNode = d.tolSameNode;
            
            DirCoarse{1}.domain    = @(coor) abs(coor(:,1) - minx) < tolSameNode;
            DirCoarse{1}.direction = [1,2];
            DirCoarse{1}.value     = 0;
            
            DirCoarse{2}.domain    = @(coor) abs(coor(:,1) - maxx) < tolSameNode;
            DirCoarse{2}.direction = [1,2];
            DirCoarse{2}.value     = 0;
            
            s.DirCond = DirCoarse;

            % disp('Nodes fixats a la mesh coarse (després Canvi 4):')
            % for i = 1:numel(s.DirCond)
            %     nodes = s.DirCond{i}.domain(s.mesh.coord);
            %     fprintf('  Dir %d: %d nodes\n', i, sum(nodes));
            % end
            % FI CANVI 4

            eifem         = EIFEM(s);
        
            % eifem.updateCoarseStiffness(alpha);
        
            ss.ddDofManager = obj.createDomainDecompositionDofManager(d);
            ss.EIFEMsolver  = eifem;
            ss.bcApplier    = bcApplierCurrent;
            ss.dMesh        = d.discMesh;
            ss.type         = 'EIFEM';
        
            eP = Preconditioner.create(ss);
            Meifem = @(r) eP.apply(r);

        end

        function mCoarse = createCoarseMesh(obj,mR,nSubdomains,tolSameNode)
            s.nsubdomains   = nSubdomains;
            s.meshReference = obj.createReferenceCoarseMesh(mR);
            s.tolSameNode   = tolSameNode;
            mRVECoarse      = MeshCreatorFromRVE.create(s);
            [mCoarse,~,~]   = mRVECoarse.create();
        end

        function cMesh = createReferenceCoarseMesh(~,mR)
            xmax = max(mR.coord(:,1));
            xmin = min(mR.coord(:,1));
            ymax = max(mR.coord(:,2));
            ymin = min(mR.coord(:,2));
        
            coord(1,1) = xmin;  coord(1,2) = ymin;
            coord(2,1) = xmax;  coord(2,2) = ymin;
            coord(3,1) = xmax;  coord(3,2) = ymax;
            coord(4,1) = xmin;  coord(4,2) = ymax;
        
            connec = [2 3 4 1];
        
            s.coord  = coord;
            s.connec = connec;

            %disp('Coords mesh coarse referència:')
            %disp(s.coord)
           
                    
            cMesh    = Mesh.create(s);
        end

        function d = createDomainDecompositionDofManager(~,data)
            s.nSubdomains             = data.nSubdomains;
            s.interfaceConnec         = data.iC;
            s.interfaceConnecReshaped = data.iCR;
            s.locGlobConnec           = data.lG;
            s.nBoundaryNodes          = data.bS{1}.mesh.nnodes;
            s.nReferenceNodes         = data.referenceMesh.nnodes;
            s.nNodes                  = data.nNodes;
            s.nDimf                   = data.nDimf;
            d = DomainDecompositionDofManager(s);
        end

        function bcApplierCurrent = createCurrentBCApplier(obj,bc)
            s.mesh               = obj.eifemData.mesh;
            s.boundaryConditions = bc;
            bcApplierCurrent     = BCApplier(s);
        end


        

    end

end