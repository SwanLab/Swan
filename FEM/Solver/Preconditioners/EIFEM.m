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
        end

        function x = apply(obj,r)
            Fcoarse = obj.projectExternalForce(r);
            x=obj.D\r;
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.mesh = cParams.mesh;
            obj.RVE  = cParams.RVE;
            obj.Kel  = repmat(obj.RVE.Kcoarse,[1,1,obj.mesh.nelem]);
            obj.DirCond = cParams.DirCond;
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
            s.fun    = []; % !!!
            trial  = LagrangianFunction.create(obj.mesh, obj.RVE.ndimf,'P1');
            test   = trial;
            assembler = AssemblerFun(s);
            LHS = assembler.assemble(obj.Kel, test, trial);
        end

        function createBoundaryConditions(obj)
            dirichlet = DirichletCondition(obj.mesh,obj.DirCond);

            s.pointloadFun = [];
            s.dirichletFun = dirichlet;
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

        function [dirichlet,pointload] = createBc(obj,boundaryMesh,dirchletBc,newmanBc)
            dirichlet = obj.createBondaryCondition(boundaryMesh,dirchletBc);
            pointload = obj.createBondaryCondition(boundaryMesh,newmanBc);
        end

        function cond = createBondaryCondition(obj,bM,condition)
            if ~isempty(condition)
                nbound = length(condition.boundaryId);
                cond = zeros(1,3);
                for ibound=1:nbound
                    ncond  = length(condition.dof(nbound,:));
                    nodeId= unique(bM{condition.boundaryId(ibound)}.globalConnec);
                    nbd   = length(nodeId);
                    for icond=1:ncond
                        bdcond= [nodeId, repmat(condition.dof(icond),[nbd,1]), repmat(condition.value(icond),[nbd,1])];
                        cond=[cond;bdcond];
                    end
                end
                cond = cond(2:end,:);
            else
                cond= [];
            end
        end

        function Fcoarse = projectExternalForce(obj,Ffine)
            nelem = obj.mesh.nelem;
            Udef  = obj.RVE.Udef;
            Urb   = obj.RVE.Ubr;
            Ut    = (Udef + Urb)';
            for ielem = nelem
                Fcoarse(ielem) = Ut*Ffine(ielem);
            end
        end

    end

end