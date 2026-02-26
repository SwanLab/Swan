classdef PreconditionerEIFEMcontinous < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        
    end
    
    properties (Access = private)
        EIFEMfilename
        weight
        nSubdomains  
        EIFEMsolver
        ddDofManager
    end
    
     properties (Access = private)
        coarseMesh
        LHS
        bcApplier        
        dir
        KeifemContinuous
        EIFEMprojection
    end    
    
    methods (Access = public)
        
         function obj = PreconditionerEIFEMcontinous(cParams)
            obj.init(cParams);
%             obj.createEIFEM();
        end
        
        
        function z = solveEIFEM(obj,r)
            RGsbd = obj.computeSubdomainResidual(r);
            uSbd =  obj.EIFEMsolver.apply(RGsbd);
            z = obj.computeContinousField(uSbd);
        end

         function K = computeKEIFEMglobal(obj,LHS)
            U = obj.computeAssembledProjectionMatrix();
            %             U = U(:,9:end);
            obj.EIFEMprojection = obj.computeReducedProjection(U);
%             LHS = obj.bcApplier.fullToReducedMatrixDirichlet(obj.LHS);
%             obj.KeifemContinuous = obj.EIFEMprojection'*LHS*obj.EIFEMprojection;
            obj.KeifemContinuous = U'*LHS*U;

        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.LHS          = cParams.LHS;
%             obj.ddDofManager = cParams.ddDofManager;
            obj.nSubdomains  = cParams.nSubdomains;
%             obj.coarseMesh   = cParams.coarseMesh;
%             obj.bcApplier     = cParams.bcApplier;
%             obj.dir = cParams.dir;
%             %             obj.EIFEMfilename = '/home/raul/Documents/Thesis/EIFEM/RAUL_rve_10_may_2024/EXAMPLE/EIFE_LIBRARY/DEF_Q4porL_2s_1.mat';
%             obj.EIFEMfilename = 'DEF_Q4porL_1.mat';
%             %           obj.EIFEMfilename = '/home/raul/Documents/Thesis/EIFEM/05_HEXAG2D/EIFE_LIBRARY/DEF_Q4auxL_1.mat';
%             obj.weight       = 0.5;
            obj.ddDofManager = cParams.ddDofManager;
            obj.EIFEMsolver  = cParams.EIFEMsolver;
            obj.bcApplier    = cParams.bcApplier;
            obj.weight       = 0.5;
%             obj.dMesh        = cParams.dMesh;
        end

%         function createEIFEM(obj)
%             filename        = obj.EIFEMfilename;
%             RVE             = TrainedRVE(filename);
%             s.RVE           = RVE;
%             s.mesh          = obj.coarseMesh;
%             s.DirCond       = obj.dir;
%             eifem           = EIFEM(s);
%             obj.EIFEMsolver = eifem;
%         end

        function uInt = computeInterfaceDisp(obj,u)
            nint = size(obj.ddDofManager.interfaceDof,3);
            uInt = zeros(size(obj.ddDofManager.interfaceDof,1),nint);
            w = [obj.weight,1- obj.weight];
            for iint = 1:nint
                ndom = size(obj.ddDofManager.interfaceDof(:,:,iint),2);
                for idom = 1:ndom
                    dom = obj.ddDofManager.interfaceDom(iint,idom);
                    unodal = u(:,dom);
                    dof = obj.ddDofManager.interfaceDof(:,idom,iint);
                    uInt(:,iint) = uInt(:,iint) + w(idom)*unodal(dof);
                end
            end
        end

        function u = smoothDisplacement(obj,u,uInterface)
            uG = obj.ddDofManager.local2global(u);
            uG = sum(uG,2);
            u  = obj.updateInterfaceValues(uG,uInterface);
        end

        function u = updateInterfaceValues(obj,u,uInterface)
            nint = size(obj.ddDofManager.interfaceDof,3);
            for iint = 1:nint
                dom = obj.ddDofManager.interfaceDom(iint,1);
                dof = obj.ddDofManager.interfaceDof(:,1,iint);
                row = ceil(dom/obj.nSubdomains(1));
                col = dom-(row-1)*obj.nSubdomains(1);
                locGlobConnec = obj.ddDofManager.localGlobalDofConnec{row,col};
                [~,ind]    = ismember(dof,locGlobConnec(:,2));
                dofGlob    = locGlobConnec(ind,1);
                u(dofGlob) = uInterface(:,iint);
            end
        end

        function R = scaleInterfaceValues(obj,R)
            nint = size(obj.ddDofManager.interfaceDof,3);
            uInt = zeros(size(obj.ddDofManager.interfaceDof,1),nint);
            w = [obj.weight,1-obj.weight];
            for iint = 1:nint
                ndom = size(obj.ddDofManager.interfaceDof(:,:,iint),2);
                for idom = 1:ndom
                    dom = obj.ddDofManager.interfaceDom(iint,idom);
                    Rnodal = R(:,dom);
                    dof = obj.ddDofManager.interfaceDof(:,idom,iint);
                    R(dof,dom) = w(idom)* R(dof,dom);
                end
            end
        end

        function Rsbd = computeSubdomainResidual(obj,R)
            RG    = obj.bcApplier.reducedToFullVectorDirichlet(R);
            %             obj.plotSolution(RG,obj.meshDomain,0,1,iter,1)
            RGsbd = obj.ddDofManager.global2local(RG);
            Rsbd = obj.scaleInterfaceValues(RGsbd);
        end

        function fc = computeContinousField(obj,f)
            fInt  = obj.computeInterfaceDisp(f);
            fc    = obj.smoothDisplacement(f,fInt);
            fc    = obj.bcApplier.fullToReducedVectorDirichlet(fc);
        end

        function Ug = computeAssembledProjectionMatrix(obj)
            Ueifem = obj.EIFEMsolver.getProjectionMatrix();
            Ueifem = repmat(Ueifem,1,obj.nSubdomains(1)*obj.nSubdomains(2));
            Ueifem = obj.scaleProjectionMatrix(Ueifem);
            Ug = obj.assembleProjectionMatrix(Ueifem);
        end

        function U = scaleProjectionMatrix(obj,U)
            nDofcoarse = size(U,2)/(obj.nSubdomains(1)*obj.nSubdomains(2));
            for i=1: nDofcoarse
                Usbd = U(:,i:nDofcoarse:end);
%                 Usbd = obj.scaleInterfaceValues(Usbd);
                Usbd = obj.ddDofManager.scaleInterfaceValues(Usbd,obj.weight);
                U(:,i:nDofcoarse:end)= Usbd;
            end
        end


        function Ug = assembleProjectionMatrix(obj,U)
            nsbd = obj.nSubdomains(1)*obj.nSubdomains(2);
            nDofcoarse = size(U,2)/nsbd;
            ndof = size(obj.LHS,1);
            Ug = zeros(ndof,nDofcoarse*nsbd);
            for i=1: nDofcoarse
                Usbd = U(:,i:nDofcoarse:end);
                Uaux = obj.ddDofManager.local2global(Usbd);
                Ug(:,i:nDofcoarse:end)= Uaux;
            end
            Ug = sparse(Ug);
        end

%         function computeKEIFEMglobal(obj)
%             U = obj.computeAssembledProjectionMatrix();
%             %             U = U(:,9:end);
%             obj.EIFEMprojection = obj.computeReducedProjection(U);
%             LHS = obj.bcApplier.fullToReducedMatrixDirichlet(obj.LHS);
%             obj.KeifemContinuous = obj.EIFEMprojection'*LHS*obj.EIFEMprojection;
% 
%         end

        function Ured = computeReducedProjection(obj,U)
            ncol = size(U,2);
            for i=1:ncol
                Ured(:,i) = obj.bcApplier.fullToReducedVectorDirichlet(U(:,i));
            end
        end

        function z = solveMEIFEMcontinuous(obj,r)
            lhs=obj.KeifemContinuous;
            phi=obj.EIFEMprojection;
            %              phi = obj.bcApplier.fullToReducedMatrixDirichlet(phi);
            r1=phi'*r;
            zP=lhs\r1;
            z =phi*zP;
            %               z = (r-LHS*z);

            %M = obj.MmodalPrecond;
            %z = M*r;
        end

        
    end
    
end