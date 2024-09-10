classdef EIFEM < handle

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
        iter
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
            obj.iter=1;
        end

        function u = apply(obj,r)
            Fcoarse = obj.projectExternalForce(r);            
            RHS     = obj.assembleRHSvector(Fcoarse);
%             RHS     = RHS + obj.reactions;
            LHSred = obj.bcApplier.fullToReducedMatrixDirichlet(obj.LHS);
            RHSred = obj.bcApplier.fullToReducedVectorDirichlet(RHS);
%             RHSred(1:2:end) = 0;
%             RHSred(2:2:end) = -0.2;
            uRed = LHSred\RHSred;
            uCoarse = obj.bcApplier.reducedToFullVectorDirichlet(uRed);
%             obj.plotSolution(uCoarse,obj.mesh,100,1,obj.iter,0)
            u = obj.reconstruct3Dsolution(uCoarse);
            obj.iter=obj.iter+1;
        end

        function u = applySubdomainNeumannDeformational(obj,r)
            Fmodal = obj.RVE.PhiDef'*r;
            uModal = obj.Kmodal\Fmodal;
            u      = obj.RVE.PhiDef*uModal;
        end

        function u = applySubdomainNeumannRigidBody(obj,r)
            Frb = obj.RVE.PhiRb'*r;
            uRb = obj.RVE.Grb\Frb;
            u      = obj.RVE.PhiRb*uRb;
        end

         function U = getProjectionMatrix(obj)
            Udef    = obj.RVE.Udef;
            Urb     = obj.RVE.Urb;
            U      = (Udef + Urb);
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.mesh = cParams.mesh;
            obj.RVE  = cParams.RVE;
            obj.Kel  = repmat(obj.RVE.Kcoarse,[1,1,obj.mesh.nelem]);
            obj.DirCond = cParams.DirCond;
            obj.dispFun = LagrangianFunction.create(obj.mesh, obj.RVE.ndimf,'P1');
            Kfine  = cParams.Kfine;
            obj.Kmodal = obj.RVE.PhiDef'*Kfine*obj.RVE.PhiDef;
        end

        function dim = getDims(obj)
            d.ndimf     = obj.RVE.ndimf;
            d.nnodes    = size(obj.mesh.coord, 1);
            d.ndofs     = d.nnodes*d.ndimf;
            d.nnodeElem = obj.mesh.nnodeElem; % should come from interp..
            d.ndofsElem = d.nnodeElem*d.ndimf;
            dim = d;
        end

        function LHS = computeLHS(obj)
            LHS = obj.assembleMatrix();
%             s.dim          = obj.getDims();
%             s.nnodeEl      = obj.mesh.nnodeElem;
%             s.globalConnec = obj.mesh.connec;
%             assembler = Assembler(s);
%             LHS = assembler.assemble(obj.Kel);
        end

        function LHS = assembleMatrix(obj)
            s.fun  = obj.dispFun; % !!!
            trial  = s.fun;
            test   = trial;
            obj.assembler = AssemblerFun(s);
            LHS = obj.assembler.assemble(obj.Kel, test, trial);
        end

        function RHS = assembleRHSvector(obj,F)
            Fcoarse = reshape(F,1,obj.dispFun.nDofsElem,[]);
            Fcoarse = permute(Fcoarse,[2 1 3]);
            fun     = [];
            RHS     = obj.assembler.assembleV(Fcoarse,fun);

%             intConec = reshape(obj.interfaceConnec',2,obj.interfaceMeshReference{1}.mesh.nnodes,[]);
%             intConec = permute(intConec,[2 1 3]);
%             Fcoarse = permute(reshape(Fcoarse',));
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
%             dirichletFun = DirichletCondition(obj.mesh, obj.DirCond{1});

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

%         function [dirichlet,pointload] = createBc(obj,boundaryMesh,dirchletBc,newmanBc)
%             dirichlet = obj.createBondaryCondition(boundaryMesh,dirchletBc);
%             pointload = obj.createBondaryCondition(boundaryMesh,newmanBc);
%         end
% 
%         function cond = createBondaryCondition(obj,bM,condition)
%             if ~isempty(condition)
%                 nbound = length(condition.boundaryId);
%                 cond = zeros(1,3);
%                 for ibound=1:nbound
%                     ncond  = length(condition.dof(nbound,:));
%                     nodeId= unique(bM{condition.boundaryId(ibound)}.globalConnec);
%                     nbd   = length(nodeId);
%                     for icond=1:ncond
%                         bdcond= [nodeId, repmat(condition.dof(icond),[nbd,1]), repmat(condition.value(icond),[nbd,1])];
%                         cond=[cond;bdcond];
%                     end
%                 end
%                 cond = cond(2:end,:);
%             else
%                 cond= [];
%             end
%         end

%         function 
%         end

        function Fcoarse = projectExternalForce(obj,Ffine)
            Udef    = obj.RVE.Udef;
            Urb     = obj.RVE.Urb;
            Ut      = (Udef + Urb)';
            Fcoarse = Ut*Ffine;
        end

       

        function u = reconstruct3Dsolution(obj,uCoarse)
            nelem = obj.mesh.nelem;
            Udef    = obj.RVE.Udef;
            Urb     = obj.RVE.Urb;
            U      = Udef + Urb;
            dofConec = obj.dispFun.getConnec();
            for ielem = 1:nelem
                uCelem = uCoarse(dofConec(ielem,:));
                u(:,ielem) =  U*uCelem;
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