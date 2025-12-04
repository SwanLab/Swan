classdef EIFEMnonPeriodic < handle

    properties (Access = public)

    end

    properties (Access = private)
        RVE
        mesh
        DirCond
    end

    properties (Access = private)
        LHS
        Kel
        boundaryConditions
        bcApplier
        assembler
        dispFun
        reactions
        LHSintegrator
        material
        meshRef
        U
        mu
        iter
    end

    methods (Access = public)

        function obj = EIFEMnonPeriodic(cParams)
            obj.init(cParams)
            obj.LHS     = obj.computeLHS(obj.mu);
            obj.computeDownscaling(obj.mu);
%             LHS = obj.computeLHS();
%             obj.LHS = LHS;
            obj.createBoundaryConditions();
            ss.mesh                 = obj.mesh;
            ss.boundaryConditions   = obj.boundaryConditions;
            obj.bcApplier           = BCApplier(ss);
%             obj.reactions           = obj.computeReactions();
        end

        function u = apply(obj,r)
            Fcoarse = obj.projectExternalForce(r);            
            RHS     = obj.assembleRHSvector(Fcoarse);

            LHSred  = obj.bcApplier.fullToReducedMatrixDirichlet(obj.LHS);
            RHSred  = obj.bcApplier.fullToReducedVectorDirichlet(RHS);
            uRed = LHSred\RHSred;
            uCoarse = obj.bcApplier.reducedToFullVectorDirichlet(uRed);
%             obj.plotSolution(uCoarse,obj.mesh,100,1,obj.iter,0)
            u = obj.reconstructSolution(uCoarse);
%                         obj.plotSolution(u(:),obj.meshRef,5,1,obj.iter,0)
                        obj.iter = obj.iter+1;
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.mesh    = cParams.mesh;
            obj.meshRef = cParams.meshRef;
            obj.RVE     = cParams.RVE;
%             obj.Kel     = repmat(obj.RVE.Kcoarse,[1,1,obj.mesh.nelem]);
            obj.DirCond = cParams.DirCond;
%             obj.dispFun = LagrangianFunction.create(obj.mesh, obj.mesh.ndim,'P2Q8');
            obj.dispFun = LagrangianFunction.create(obj.mesh, obj.mesh.ndim,'P1');
%             obj.LHSintegrator = obj.createLHSintegrator();
            if length(cParams.mu) == 1
                obj.mu = cParams.mu*ones(cParams.mesh.nelem,1);
            else
                obj.mu = reshape(cParams.mu',1,[]);
            end
            obj.iter=1;
        end

        function computeElementalLHS(obj,mu)
            nelem = obj.mesh.nelem;
            for i = 1:nelem
                obj.Kel(:,:,i) = obj.RVE.Kcoarse(mu(i));
%                   obj.Kel(:,:,i) = obj.RVE.Kcoarse(:,:,i);
            end
        end

        function computeDownscaling(obj,mu)
            nelem = obj.mesh.nelem;
            for i = 1:nelem
                Udef = obj.RVE.Udef(mu(i));
                Urb  = obj.RVE.Urb(mu(i));
%                 Udef = obj.RVE.Udef(:,:,i);
%                 Urb  = obj.RVE.Urb(:,:,i);
                obj.U(:,:,i) = Udef + Urb;
            end
        end

        function LHSint = createLHSintegrator(obj)
            s.type     = 'ElasticStiffnessMatrixModal';
            s.mesh     = obj.mesh;
            s.test     = [];
            s.trial    = [];
            s.quadratureOrder = 2;
            LHSint = LHSIntegrator.create(s);
        end

        function T = computeDownscalingFunction(obj)
            nElem = obj.mesh.nelem;
            functionType = 'P1';
            for i = 1:nElem
             Udef    = obj.RVE.Udef;
             Urb     = obj.RVE.Urb;
             U       = Urb + Udef; 
             T = ModalFunction.create(obj.meshRef,U,functionType);
            end
        end

        function LHS = computeLHS(obj,mu)
            obj.computeElementalLHS(mu)
            LHS = obj.assembleMatrix(obj.Kel,obj.dispFun,obj.dispFun);
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
                dir = DirichletCondition(obj.mesh, obj.DirCond{i},obj.dispFun.order);
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
            nelem = obj.mesh.nelem;
            for ielem = 1:nelem
                Fcoarse(:,ielem) = obj.U(:,:,ielem)'*Ffine(:,ielem);
            end
% 
%             Udef    = obj.RVE.Udef(mu);
%             Urb     = obj.RVE.Urb(mu);
%             Ut      = (Udef + Urb)';
%             Fcoarse = Ut*Ffine;
        end

        function u = reconstructSolution(obj,uCoarse)
            nElem = obj.mesh.nelem;
%             Udef  = obj.RVE.Udef;
%             Urb   = obj.RVE.Urb;
%             U     = Udef + Urb;
            dofConec = obj.dispFun.getDofConnec();
            for ielem = 1:nElem
                uCelem = uCoarse(dofConec(ielem,:));
                u(:,ielem) =  obj.U(:,:,ielem)*uCelem;
            end
        end

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