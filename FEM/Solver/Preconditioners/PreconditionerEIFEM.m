classdef PreconditionerEIFEM < handle

    properties (Access = public)

    end

    properties (Access = private)
        EIFEMfilename
        weight
        EIFEMsolver
    end

    properties (Access = private)
        coarseMesh
        LHS
        bcApplier
        dir
        ddDofManager
        nSubdomains
    end

    methods (Access = public)

        function obj = PreconditionerEIFEM(cParams)
            obj.init(cParams);
            obj.createEIFEM();
        end

        function z = solveEIFEM(obj,r)
            RGsbd = obj.computeSubdomainResidual(r);
            uSbd =  obj.EIFEMsolver.apply(RGsbd);
            z = obj.computeContinousField(uSbd);
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.LHS          = cParams.LHS;
            obj.ddDofManager = cParams.ddDofManager;
            obj.coarseMesh   = cParams.coarseMesh;
            obj.bcApplier     = cParams.bcApplier;
            obj.dir = cParams.dir;
            % obj.EIFEMfilename = '/home/raul/Documents/Thesis/EIFEM/RAUL_rve_10_may_2024/EXAMPLE/EIFE_LIBRARY/DEF_Q4porL_2s_1.mat';
            obj.EIFEMfilename = 'DEF_Q4porL_1.mat';
            % obj.EIFEMfilename = '/home/raul/Documents/Thesis/EIFEM/05_HEXAG2D/EIFE_LIBRARY/DEF_Q4auxL_1.mat';
            obj.weight       = 0.5;
        end

        function createEIFEM(obj)
            filename        = obj.EIFEMfilename;
            RVE             = TrainedRVE(filename);
            s.RVE           = RVE;
            s.mesh          = obj.coarseMesh;
            s.DirCond       = obj.dir;
            eifem           = EIFEM(s);
            obj.EIFEMsolver = eifem;
        end

%         function uInt = computeInterfaceDisp(obj,u)
%             nint = size(obj.ddDofManager.interfaceDof,3);
%             uInt = zeros(size(obj.ddDofManager.interfaceDof,1),nint);
%             w = [obj.weight,1- obj.weight];
%             for iint = 1:nint
%                 ndom = size(obj.ddDofManager.interfaceDof(:,:,iint),2);
%                 for idom = 1:ndom
%                     dom = obj.ddDofManager.interfaceDom(iint,idom);
%                     unodal = u(:,dom);
%                     dof = obj.ddDofManager.interfaceDof(:,idom,iint);
%                     uInt(:,iint) = uInt(:,iint) + w(idom)*unodal(dof);
%                 end
%             end
%         end

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

%         function u = smoothDisplacement(obj,u,uInterface)
%             uG = obj.ddDofManager.local2global(u);
%             uG = sum(uG,2);
%             u  = obj.updateInterfaceValues(uG,uInterface);
%         end
% 
%         function u = updateInterfaceValues(obj,u,uInterface)
%             nint = size(obj.ddDofManager.interfaceDof,3);
%             nSub = obj.ddDofManager.nSubdomains(1);
%             for iint = 1:nint
%                 dom = obj.ddDofManager.interfaceDom(iint,1);
%                 dof = obj.ddDofManager.interfaceDof(:,1,iint);
%                 row = ceil(dom/nSub);
%                 col = dom-(row-1)*nSub;
%                 locGlobConnec = obj.ddDofManager.localGlobalDofConnec{row,col};
%                 [~,ind]    = ismember(dof,locGlobConnec(:,2));
%                 dofGlob    = locGlobConnec(ind,1);
%                 u(dofGlob) = uInterface(:,iint);
%             end
%         end

        function RGsbd = computeSubdomainResidual(obj,R)
            RG    = obj.bcApplier.reducedToFullVectorDirichlet(R);
            %             obj.plotSolution(RG,obj.meshDomain,0,1,iter,1)
            RGsbd = obj.ddDofManager.global2local(RG);  %dissemble
            RGsbd = obj.scaleInterfaceValues(RGsbd);    %scale
        end

        function fc = computeContinousField(obj,f)
%             fInt  = obj.computeInterfaceDisp(f);
            fS  = obj.scaleInterfaceValues(f);         %scale
            fG  = obj.ddDofManager.local2global(fS);   %assemble
            fc  = sum(fG,2);                           %assemble
%             fc    = obj.smoothDisplacement(f,fInt);
            fc  = obj.bcApplier.fullToReducedVectorDirichlet(fc);
        end

    end

end