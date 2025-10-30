classdef EIFEM < handle

    properties (Access = public)

    end

    properties (Access = private)
        RVE
        mesh
        DirCond
        weight
    end

    properties (Access = private)
        LHS
        Kel
        boundaryConditions
        ddDofManager
        bcApplier
        assembler
        dispFun
        Kmodal
        reactions
    end

    methods (Access = public)

        function obj = EIFEM(cParams)
            obj.init(cParams)
            LHS = obj.computeLHS();
            obj.LHS = LHS;
            obj.createBoundaryConditions();
            ss.mesh                 = obj.mesh;
            ss.boundaryConditions   = obj.boundaryConditions;
            obj.bcApplier           = BCApplier(ss);
            obj.reactions           = obj.computeReactions();
        end

        function u = apply(obj,r)
            Fcoarse = obj.projectExternalForce(r);            
            RHS     = obj.assembleRHSvector(Fcoarse);
            LHSred = obj.bcApplier.fullToReducedMatrixDirichlet(obj.LHS);
            RHSred = obj.bcApplier.fullToReducedVectorDirichlet(RHS);
            uRed = LHSred\RHSred;
%             uCoarse = obj.bcApplier.reducedToFullVectorDirichlet(uRed);
%             obj.plotSolution(uCoarse,obj.mesh,100,1,obj.iter,0)
            u = obj.reconstructSolution(uRed);
        end

        function M = reduceMatrix(obj,M)
           M =  obj.bcApplier.fullToReducedMatrixDirichlet(M);
        end       

        function u = reconstructSolution(obj,uCoarse)
            uCoarse = obj.bcApplier.reducedToFullVectorDirichlet(uCoarse);
            nElem = obj.mesh.nelem;
            Udef  = obj.RVE.Udef;
            Urb   = obj.RVE.Urb;
            U     = Udef + Urb;
            dofConec = obj.dispFun.getDofConnec();
            for ielem = 1:nElem
                uCelem = uCoarse(dofConec(ielem,:));
                u(:,ielem) =  U*uCelem;
            end
        end
           
        function A = assembleMatrix(obj,Aelem,f1,f2)
            dofsF1 = f1.getDofConnec();
            if isequal(f1, f2)
                dofsF2 = dofsF1;
            else
                dofsF2 = f2.getDofConnec();
            end
            nDofs1     = numel(f1.fValues);
            nDofs2     = numel(f2.fValues);
            ndofsElem1 = size(Aelem, 1);
            ndofsElem2 = size(Aelem, 2);

            [iElem, jElem] = meshgrid(1:ndofsElem1, 1:ndofsElem2);
            iElem = iElem(:);
            jElem = jElem(:);

            dofsI = dofsF1(:, iElem);
            dofsJ = dofsF2(:, jElem);

            rowIdx = dofsI(:);
            colIdx = dofsJ(:);
            Aval   = permute(Aelem,[3 2 1]);
            values = Aval(:);
            A = sparse(rowIdx, colIdx, values, nDofs1, nDofs2);
        end

        function uC = computeContinousField(obj,uD)
            fS  = obj.ddDofManager.scaleInterfaceValues(uD,obj.weight);         %scale
            fG  = obj.ddDofManager.local2global(fS);   %assemble
            uC  = sum(fG,2);                           %assemble
            %uC  = obj.bcApplier.fullToReducedVectorDirichlet(uC);
        end

        function uFine = injectCoarseToFineBoundary(obj, uCoarse)
            % 1) undo Dirichlet reduction on the *coarse* vector
            uCfull = obj.bcApplier.reducedToFullVectorDirichlet(uCoarse);
        
            % 2) build a local (discontinuous) vector with ONLY interface/boundary dofs set
            %    (same shape you pass to ddDofManager.local2global later)
            %    Start by expanding the *global* coarse values to local-by-subdomain copies:
            %uD = obj.ddDofManager.global2local(uCfull);         % scatter to subdomains
            uD = obj.ddDofManager.scaleInterfaceValues(uD, 0.5);% optional: same weights you use elsewhere
        
            % 3) assemble the discontinuous copies back to the global fine vector
            fG  = obj.ddDofManager.local2global(uD);            % assemble (sum into global)
            uFine = sum(fG,2);                                   % global fine vector length = discMesh.ndim*discMesh.nnodes
        
            % 4) zero out interior fine DOFs if needed (optional):
            %    keep only boundary nodes; if you have a mask of interface dofs, apply it here
            % mask = obj.ddDofManager.interfaceMaskGlobal();     % (illustrative)
            % uFine(~mask) = 0;
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.mesh    = cParams.mesh;
            obj.RVE     = cParams.RVE;
            obj.Kel     = repmat(obj.RVE.Kcoarse,[1,1,obj.mesh.nelem]);
            obj.DirCond = cParams.DirCond;
%             obj.dispFun = LagrangianFunction.create(obj.mesh, obj.RVE.ndimf,'P1');
            obj.dispFun = LagrangianFunction.create(obj.mesh, obj.mesh.ndim,'P1');
            obj.ddDofManager = cParams.ddDofManager;
            obj.weight       = 0.5;
        end

        function LHS = computeLHS(obj)
            LHS = obj.assembleMatrix(obj.Kel,obj.dispFun,obj.dispFun);
        end
        
      

        function RHS = assembleRHSvector(obj,F)
            Fcoarse = reshape(F,1,obj.dispFun.nDofsElem,[]);
            Fcoarse = squeeze(Fcoarse);
            RHS     = obj.assembleVector(Fcoarse,obj.dispFun);
        end

        function F = assembleVector(obj,Felem, f)
            dofConnec = f.getDofConnec();
            nDofs     = numel(f.fValues);
            rowIdx    = dofConnec(:);
            Felem = Felem';
            F = sparse(rowIdx, 1, Felem(:), nDofs, 1);
        end

        function R = computeReactions(obj)
            bc      = obj.boundaryConditions;
            dirich  = bc.dirichlet_dofs;
            dirichV = bc.dirichlet_vals;
            if ~isempty(dirich)
                R = -obj.LHS(:,dirich)*dirichV;
            else
                R = zeros(sum(obj.dim.ndofs(:)),1);
            end
        end

        function createBoundaryConditions(obj)
            dirichletFun = [];
             for i = 1:numel(obj.DirCond)
                dir = DirichletCondition(obj.mesh, obj.DirCond{i});
                dirichletFun = [dirichletFun, dir];
            end

            s.pointloadFun = [];
            s.dirichletFun = dirichletFun;
            s.periodicFun  = [];
            s.mesh         = obj.mesh;
            bc             = BoundaryConditions(s);
%             bM = obj.mesh.createBoundaryMesh();

%             dirichletBc.boundaryId=1;
%             dirichletBc.dof=[1,2];
%             dirichletBc.value=[0,0];
%             %             newmanBc.boundaryId=2;
%             %             newmanBc.dof=[2];
%             %             newmanBc.value=[10];
%             newmanBc= [];
% 
%             [dirichlet,pointload] = obj.createBc(bM,dirichletBc,newmanBc);
%             bc.dirichlet=dirichlet;
%             bc.pointload=pointload;
            obj.boundaryConditions = bc;
        end

        function Fcoarse = projectExternalForce(obj,Ffine)
            Udef    = obj.RVE.Udef;
            Urb     = obj.RVE.Urb;
            Ut      = (Udef + Urb)';
            Fcoarse = Ut*Ffine;
        end

%         function u = reconstructSolution(obj,uCoarse)
%             nElem = obj.mesh.nelem;
%             Udef  = obj.RVE.Udef;
%             Urb   = obj.RVE.Urb;
%             U     = Udef + Urb;
%             dofConec = obj.dispFun.getDofConnec();
%             for ielem = 1:nElem
%                 uCelem = uCoarse(dofConec(ielem,:));
%                 u(:,ielem) =  U*uCelem;
%             end
%         end

         function plotSolution(obj,x,mesh,row,col,iter,flag)
            if nargin <7
                 flag =0;
            end
            %             xFull = bc.reducedToFullVector(x);
            if size(x,2)==1
                 s.fValues = reshape(x,2,[])';
            else
                 s.fValues = x; 
            end
%            
            
            s.mesh = mesh;
            s.fValues(:,end+1) = 0;
            s.ndimf = 2;
            s.order = obj.dispFun.order;
            xF = LagrangianFunction(s);
%             xF.plot();
            if flag == 0
                xF.print(['domain',num2str(row),num2str(col),'_',num2str(iter)],'Paraview')
            elseif flag == 1
                xF.print(['DomainResidual',num2str(row),num2str(col),'_',num2str(iter)],'Paraview')
            elseif flag == 2
                xF.print(['Residual',num2str(row),num2str(col),'_',num2str(iter)],'Paraview')
            elseif flag == 3
                xF.print(['domainFine',num2str(row),num2str(col),'_',num2str(iter)],'Paraview')
            elseif flag == 4
                xF.print(['domainNeuman',num2str(row),num2str(col),'_',num2str(iter)],'Paraview')
            end
            fclose('all');
        end


    end

end