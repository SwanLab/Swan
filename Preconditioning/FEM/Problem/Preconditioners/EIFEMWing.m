classdef EIFEMWing < handle

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
        Kmodal
        reactions
    end

    methods (Access = public)

        function obj = EIFEMWing(cParams)
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
            uCoarse = obj.bcApplier.reducedToFullVectorDirichlet(uRed);
%             obj.plotSolution(uCoarse,obj.mesh,100,1,obj.iter,0)
            u = obj.reconstructSolution(uCoarse);
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
        end

        function LHS = computeLHS(obj)
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
            F = sparse(Felem);
        end

        function R = computeReactions(obj)
            bc      = obj.boundaryConditions;
            dirich  = bc.dirichlet_dofs;
            dirichV = bc.dirichlet_vals;
            if ~isempty(dirich)
                R = -obj.LHS(:,dirich)*dirichV;
            else
                % Usar nDofs del dispFun en lugar de obj.dim.ndofs
                nDofs = obj.dispFun.nDofs;
                R = zeros(nDofs, 1);
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
            
            nU_rows = size(obj.RVE.Udef, 1);
            max_dof = max([obj.RVE.DOFl(:); obj.RVE.DOFr(:)]);
            nDOF_fine = length(obj.RVE.DOFl) + length(obj.RVE.DOFr);
            % r está en espacio fine reducido (DOFl)
            % 1. Expandir a espacio fine completo (asumiendo DOFr = 0 para residuo)
            % El residuo en DOFs restringidos es cero por definición
            
            % Crear r_full con el tamaño correcto
            if nU_rows >= max_dof
                % U tiene suficiente espacio (puede incluir DOFs adicionales)
                r_full = zeros(nU_rows, size(Ffine, 2));
                % Mapear DOFl a las posiciones correctas
                % DOFl son índices que apuntan a posiciones en el espacio fine completo
                valid_indices = obj.RVE.DOFl <= nU_rows;
                if any(~valid_indices)
                    if options.verbose
                        warning('Algunos DOFl exceden dimensiones de U. Ajustando.');
                    end
                    DOFl_valid = obj.RVE.DOFl(valid_indices);
                    r_full(DOFl_valid, :) = Ffine(valid_indices, :);
                else
                    r_full(obj.RVE.DOFl, :) = Ffine;
                end
            else
                % U es más pequeño, truncar si es necesario
                r_full = zeros(nDOF_fine, size(Ffine, 2));
                r_full(obj.RVE.DOFl, :) = Ffine;
                if nU_rows < nDOF_fine
                    r_full = r_full(1:nU_rows, :);
                end
            end

            Udef    = obj.RVE.Udef;
            Urb     = obj.RVE.Urb;
            Ut      = (Udef + Urb)';
            Fcoarse = Ut*r_full;
        end

        function u = reconstructSolution(obj,uCoarse)
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