classdef Coarse < handle

    properties (Access = public)

    end

    properties (Access = private)
        RVE
        mesh
        DirCond
    end

    properties (Access = private)
        LHS
        U
        Kel
        boundaryConditions
        bcApplier
        assembler
        dispFun
        Kmodal
        reactions
    end

    methods (Access = public)

        function obj = Coarse(cParams)
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
            RHS     = obj.assembleRHSvector1(Fcoarse);
            LHSred = obj.bcApplier.fullToReducedMatrixDirichlet(obj.LHS);
            RHSred = obj.bcApplier.fullToReducedVectorDirichlet(RHS);
            uRed = LHSred\RHSred;
            uCoarse = obj.bcApplier.reducedToFullVectorDirichlet(uRed);
%             obj.plotSolution(uCoarse,obj.mesh,100,1,obj.iter,0)
            u = obj.reconstructSolution(uCoarse);
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.mesh    = cParams.mesh;
            obj.RVE     = cParams.RVE;
            %obj.Kel     = repmat(obj.RVE.Kcoarse,[1,1,obj.mesh.nelem]);
            obj.DirCond = cParams.DirCond;
            
            obj.Kel = [];
            obj.U = [];

            for i=1:size(obj.RVE,1)
                for j=1:size(obj.RVE,2)
                    obj.Kel = cat(3, obj.Kel, obj.RVE{i,j}.Kcoarse);
                    obj.U   = cat(3, obj.U, obj.RVE{i,j}.U);
                end
            end

            
            obj.dispFun = LagrangianFunction.create(obj.mesh,obj.RVE{1,1}.ndimf,'P1');
            %obj.dispFun = LagrangianFunction.create(obj.mesh, obj.RVE{1,1}.ndimf,'P1');
            
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
            %LHS = obj.assembleMatrix1();
            LHS = obj.assembleMatrix1(obj.Kel,obj.dispFun,obj.dispFun);

        end

        function LHS = assembleMatrix(obj)
            s.fun  = obj.dispFun; % !!!
            trial  = s.fun;
            test   = trial;
            obj.assembler = AssemblerFun(s);
            LHS = obj.assembler.assemble(obj.Kel, test, trial);
        end


        function A = assembleMatrix1(obj,Aelem,f1,f2)
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

        function RHS = assembleRHSvector1(obj,F)
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


        function RHS = assembleRHSvector(obj,F)
            Fcoarse = reshape(F,1,obj.dispFun.nDofsElem,[]);
            Fcoarse = permute(Fcoarse,[2 1 3]);
            fun     = [];
            RHS     = obj.assembler.assembleV(Fcoarse,fun);

        end

        function R = computeReactions(obj)
            bc      = obj.boundaryConditions;
            dirich  = bc.dirichlet_dofs;
            dirichV = bc.dirichlet_vals;
            if ~isempty(dirich)
                R = -obj.LHS(:,dirich)*dirichV;  %% 
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
           %--------------------------------------------
           Fcoarse = [];
           a = 1;
           for i = 1:size(obj.RVE,1)
                for j = 1:size(obj.RVE,2)
                    Fcoarse = cat(3, Fcoarse, obj.U(:,:,a)'*Ffine(:,a));
                    a = a+1;
                end

           end
            %Ut      = obj.U';
            %Fcoarse = Ut*Ffine;
        end

        function u = reconstructSolution(obj,uCoarse)
            nElem = obj.mesh.nelem;
            dofConec = obj.dispFun.getDofConnec();
            % for ielem = 1:nElem
            %     uCelem = uCoarse(dofConec(ielem,:));
            %     u(:,ielem) =  obj.RVE.U*uCelem;
            % end
            
            ielem = 1;
            for i = 1:size(obj.RVE,1)
                for j = 1:size(obj.RVE,2)
                    uCelem = uCoarse(dofConec(ielem,:));
                    u(:,ielem) =  obj.RVE{i,j}.U*uCelem;
                    ielem = ielem+1;
                end
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