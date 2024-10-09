classdef GeneralPreconditioner < handle

    properties (Access = public)

    end

    properties (Access = private)
        EIFEMfilename
        
     

        localGlobalDofConnec
        interfaceDof
        interfaceDom
        weight

        nSubdomains

        scale
        ndimf
        functionType
        solverCase
        
        displacementFun
        EIFEMsolver
       
    end

    properties (Access = private)
        meshDomain
        coarseMesh
        LHS
        RHS
        bcApplier
        interfaceConnec
        locGlobConnec        
        nBoundaryNodes
        nReferenceNodes
        dir
    end    

    methods (Access = public)

        function obj = GeneralPreconditioner(cParams)
            obj.init(cParams);
            obj.constructionObject()
        end

        function z = solveEIFEM(obj,r)
            RGsbd = obj.computeSubdomainResidual(r);
            uSbd =  obj.EIFEMsolver.apply(RGsbd);
            z = obj.computeContinousField(uSbd);
        end
        
        function x = InexactCG(obj,r,A,P)
            x0 = zeros(size(r));
            factor = 0.1;
            tol = factor*norm(r);
            x = PCG.solve(A,r,x0,P,tol);
            tau = 1;
            
            %tau = @(r,A) 1;
         %   tau = @(r,A) r'*r/(r'*A(r));
         %   x = RichardsonSolver.solve(A,r,x0,P,tol,tau);
            
        end
        
 
    end

    methods (Access = private)

        function init(obj,cParams)
            obj.LHS = cParams.LHS;
            obj.RHS = cParams.RHS;
            obj.nSubdomains = cParams.nSubdomains;
            obj.meshDomain  = cParams.meshDomain;
            obj.coarseMesh = cParams.coarseMesh;
            obj.interfaceConnec = cParams.interfaceConnec;
            obj.locGlobConnec = cParams.locGlobConnec;
            obj.nBoundaryNodes = cParams.nBoundaryNodes;
            obj.nReferenceNodes = cParams.nReferenceNodes;
            obj.bcApplier = cParams.bcApplier;
            obj.dir = cParams.dir;
            obj.displacementFun = LagrangianFunction.create(obj.meshDomain,obj.meshDomain.ndim,'P1');            
            obj.scale        = 'MACRO';
            obj.ndimf        = 2;
            obj.functionType = 'P1';
            obj.solverCase   = 'REDUCED';
            %             obj.EIFEMfilename = '/home/raul/Documents/Thesis/EIFEM/RAUL_rve_10_may_2024/EXAMPLE/EIFE_LIBRARY/DEF_Q4porL_2s_1.mat';
            obj.EIFEMfilename = 'DEF_Q4porL_1.mat';
            %           obj.EIFEMfilename = '/home/raul/Documents/Thesis/EIFEM/05_HEXAG2D/EIFE_LIBRARY/DEF_Q4auxL_1.mat';
            obj.weight       = 0.5;
        end

        function constructionObject(obj)

          
            obj.localGlobalDofConnec = obj.createlocalGlobalDofConnec();
            [obj.interfaceDof,obj.interfaceDom] = obj.computeLocalInterfaceDof();

          

            meshDomainCoarse = obj.coarseMesh;

            obj.EIFEMsolver = obj.createEIFEM(meshDomainCoarse,obj.dir);
            %             obj.computeKEIFEMglobal();

            

        end

         function dim = getFunDims(obj)
            d.ndimf  = obj.displacementFun.ndimf;
            d.nnodes = size(obj.displacementFun.fValues, 1);
            d.ndofs  = d.nnodes*d.ndimf;
            d.nnodeElem = obj.meshDomain.nnodeElem; % should come from interp..
            d.ndofsElem = d.nnodeElem*d.ndimf;
            dim = d;
        end

        function  localGlobalDofConnec = createlocalGlobalDofConnec(obj)
            ndimf = obj.displacementFun.ndimf;
            ndom  = obj.nSubdomains(1)*obj.nSubdomains(2);
            for dom = 1:ndom
                row = ceil(dom/obj.nSubdomains(1));
                col = dom-(row-1)*obj.nSubdomains(1);
                localGlobalDofConnecDom = zeros(1,2);
                nodeG = obj.locGlobConnec(:,1,dom);
                nodeL = obj.locGlobConnec(:,2,dom);
                for iunkn = 1:ndimf
                    dofConec = [ndimf*(nodeG - 1) + iunkn ,  ndimf*(nodeL - 1) + iunkn] ;
                    localGlobalDofConnecDom = [localGlobalDofConnecDom;dofConec ];
                end
                localGlobalDofConnec{row,col} = localGlobalDofConnecDom(2:end,:);
            end
        end

        function [interfaceDof,interfaceDom] = computeLocalInterfaceDof(obj)
            intConec = reshape(obj.interfaceConnec',2,obj.nBoundaryNodes,[]);
            intConec = permute(intConec,[2 1 3]);
            nint = size(intConec,3);
            globaldof=0;
            ndimf = obj.ndimf;
            ndofs = obj.nReferenceNodes*ndimf;

            for iint=1:nint
                ndom = size(intConec,2); %length(intConec(1,:,iint));
                for idom = 1:ndom
                    dofaux=0;
                    nodesI = intConec(:,idom,iint);
                    dom = ceil(intConec(1,idom,iint)/obj.nReferenceNodes);
                    globaldof = (dom-1)*ndofs;
                    for iunkn=1:ndimf
                        DOF = ndimf*(nodesI - 1) + iunkn;
                        DOF = DOF-globaldof;
                        dofaux= [dofaux; DOF];
                    end
                    interfaceDof(:,idom,iint) = dofaux(2:end);
                    interfaceDom(iint,idom) = dom;
                    %                     globaldof = globaldof + (iint*(idom-1)+iint)*dim.ndofs;
                end
                %                 interfaceDof(:,iint) = dofaux(2:end);
                %                 globaldof = globaldof + dim.ndofs;
            end
        end

        function Gvec = local2global(obj,Lvec)
            %             ndimf  = obj.displacementFun.ndimf;
            Gvec   = zeros(obj.displacementFun.nDofs,obj.nSubdomains(1)*obj.nSubdomains(2));
            %             Gvec(locGlobConnec(:,1)) = Lvec(locGlobConnec(:,2));
            ind    = 1;
            for jdom = 1: obj.nSubdomains(2)
                for idom = 1: obj.nSubdomains(1)
                    locGlobConnec = obj.localGlobalDofConnec{jdom,idom};
                    Gvec(locGlobConnec(:,1),ind) = Lvec(locGlobConnec(:,2),ind);
                    ind=ind+1;
                end
            end
        end

        function Lvec = global2local(obj,Gvec)
            ndimf  = obj.displacementFun.ndimf;
            Lvec   = zeros(obj.nReferenceNodes*ndimf,obj.nSubdomains(1)*obj.nSubdomains(2));
            ind    = 1;
            for jdom = 1: obj.nSubdomains(2)
                for idom = 1: obj.nSubdomains(1)
                    locGlobConnec = obj.localGlobalDofConnec{jdom,idom};
                    Lvec(locGlobConnec(:,2),ind) = Gvec(locGlobConnec(:,1));
                    ind=ind+1;
                end
            end
        end

        function eifem = createEIFEM(obj,meshDomainCoarse,Dir)
            filename   = 'DEF_Q4porL_1.mat';
            RVE        = TrainedRVE(filename);
            s.RVE      = RVE;
            s.mesh     = meshDomainCoarse;
            s.DirCond  = Dir;
            eifem      = EIFEM(s);
        end

        function uInt = computeInterfaceDisp(obj,u)
            nint = size(obj.interfaceDof,3);
            uInt = zeros(size(obj.interfaceDof,1),nint);
            w = [obj.weight,1- obj.weight];
            for iint = 1:nint
                ndom = size(obj.interfaceDof(:,:,iint),2);
                for idom = 1:ndom
                    dom = obj.interfaceDom(iint,idom);
                    unodal = u(:,dom);
                    dof = obj.interfaceDof(:,idom,iint);
                    uInt(:,iint) = uInt(:,iint) + w(idom)*unodal(dof);
                end
            end
        end

        function u = smoothDisplacement(obj,u,uInterface)
            uG = obj.local2global(u);
            uG = sum(uG,2);
            u  = obj.updateInterfaceValues(uG,uInterface);
        end

        function u = updateInterfaceValues(obj,u,uInterface)
            nint = size(obj.interfaceDof,3);
            for iint = 1:nint
                dom = obj.interfaceDom(iint,1);
                dof = obj.interfaceDof(:,1,iint);
                row = ceil(dom/obj.nSubdomains(1));
                col = dom-(row-1)*obj.nSubdomains(1);
                locGlobConnec = obj.localGlobalDofConnec{row,col};
                [~,ind]    = ismember(dof,locGlobConnec(:,2));
                dofGlob    = locGlobConnec(ind,1);
                u(dofGlob) = uInterface(:,iint);
            end
        end

        function R = scaleInterfaceValues(obj,R)
            nint = size(obj.interfaceDof,3);
            uInt = zeros(size(obj.interfaceDof,1),nint);
            w = [obj.weight,1-obj.weight];
            for iint = 1:nint
                ndom = size(obj.interfaceDof(:,:,iint),2);
                for idom = 1:ndom
                    dom = obj.interfaceDom(iint,idom);
                    Rnodal = R(:,dom);
                    dof = obj.interfaceDof(:,idom,iint);
                    R(dof,dom) = w(idom)* R(dof,dom);
                end
            end
        end

        function Rsbd = computeSubdomainResidual(obj,R)
            RG    = obj.bcApplier.reducedToFullVectorDirichlet(R);
            %             obj.plotSolution(RG,obj.meshDomain,0,1,iter,1)
            RGsbd = obj.global2local(RG);
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
            nDofcoarse = size(U,2)/obj.nSubdomains(1)*obj.nSubdomains(2);
            for i=1: nDofcoarse
                Usbd = U(:,i:nDofcoarse:end);
                Usbd = obj.scaleInterfaceValues(Usbd);
                U(:,i:nDofcoarse:end)= Usbd;
            end
        end


        function Ug = assembleProjectionMatrix(obj,U)
            nsbd = obj.nSubdomains(1)*obj.nSubdomains(2);
            nDofcoarse = size(U,2)/nsbd;
            Ug = zeros(obj.displacementFun.nDofs,nDofcoarse*nsbd);
            for i=1: nDofcoarse
                Usbd = U(:,i:nDofcoarse:end);
                Uaux = obj.local2global(Usbd);
                Ug(:,i:nDofcoarse:end)= Uaux;
            end
            Ug = sparse(Ug);
        end

        function computeKEIFEMglobal(obj)
            U = obj.computeAssembledProjectionMatrix();
            %             U = U(:,9:end);
            obj.EIFEMprojection = obj.computeReducedProjection(U);
            LHS = obj.bcApplier.fullToReducedMatrixDirichlet(obj.LHS);
            obj.KeifemContinuous = obj.EIFEMprojection'*LHS*obj.EIFEMprojection;

        end

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
            s.order = obj.functionType;
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

    methods (Access = public, Static)

        function z = multiplePrec(r,P1,P2,P3,A)
            z1 = P1(r);
            r  = r-A(z1);
            z2 = P2(r);
            r  = r-A(z2);
            z3 = P3(r);
            z  = z1+z2+z3;
        end

        function z = multiplePrec2(r,P1,P2,A)
            z1 = P1(r);
            r  = r-A(z1);
            z2 = P2(r);
            z  = z1+z2;
        end        

        function z = additivePrec(r,P1,P2)
            z1 = P1(r);
            z2 = P2(r);
            z  = z1+z2;
        end           

    end

end