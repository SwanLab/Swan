classdef BalancingNeumannNeumann_ndom < handle

    properties (Access = public)

    end

    properties (Access = private)

    end

    properties (Access = private)
        meshDomain
        boundaryConditions
        material
        meshReference
        interfaceMeshReference
        meshSubDomain
        interfaceConnec
        nSubdomains
        ninterfaces
        interfaceMeshSubDomain
        globalMeshConnecSubDomain
        interfaceMeshConnecSubDomain
        subDomainContact
        cornerNodes
        quad
        interfaceDof
        interfaceDom
        dofInterfaceDomain
        interiorDof
        hExt
        locGlobConnec
        localGlobalDofConnec
        coarseSpace
        domainFun
        domLHS

        displacementFun
        LHS
        RHS
        scale
        theta
        Fext
        weight
        interfaceNeumanDof
        interfaceDirichletDof
    end

    methods (Access = public)

        function obj = BalancingNeumannNeumann_ndom()
            close all
            obj.init();
            obj.createReferenceMesh();

            s.nsubdomains = obj.nSubdomains; %nx ny
            s.meshReference = obj.meshReference;
            m = MeshCreatorFromRVE(s);
            [obj.meshDomain,obj.meshSubDomain,obj.interfaceConnec,~,obj.locGlobConnec] = m.create();

            obj.domainFun = P1Function.create(obj.meshDomain, obj.meshDomain.ndim);
            obj.createSubdomainDisplacementFun();
            [obj.interfaceDof,obj.interfaceDom] = obj.computeLocalInterfaceDof();
            obj.dofInterfaceDomain = obj.assingInterfaceDof2Domain();
            obj.interiorDof        = obj.assingInteriorDof();
            obj.computeSubdomainLHS();
            obj.createSubdomainBoundaryConditions();

            obj.localGlobalDofConnec = obj.createlocalGlobalDofConnec();

            obj.computeDomainLHS();

            coarseSpace = obj.computeCoarseSpace();
            obj.plotfields(coarseSpace)


%             obj.createSubdomainBoundaryConditions();

            obj.interfaceNeumanDof = obj.identifyNeumanInterfaceDof();
            obj.interfaceDirichletDof = obj.identifyDirichletInterfaceDof();

            obj.NeumannNeumannSolver();


            %             obj.solveDomainProblem();



        end

        function solveDomainProblem(obj)
            s.mesh     = obj.meshDomain;
            s.bc       = obj.boundaryConditions;
            s.material = obj.material;
            s.type     = 'ELASTIC';
            s.scale    = 'MACRO';
            s.dim      = '2D';
            s.solverTyp = 'ITERATIVE';
            s.iterativeSolverTyp = 'PCG';
            s.preconditionerType = 'EIFEM';
            s.tol = 1e-6;

            fem        = FEM.create(s);
            fem.solve();
        end

    end

    methods (Access = private)

        function init(obj)
            obj.nSubdomains = [3 2]; %nx ny
            obj.scale    = 'MACRO';
            obj.weight   = 0.5;
            obj.theta    = 0.1;
        end

        function createReferenceMesh(obj)
            filename   = 'lattice_ex1';
            a.fileName = filename;
            femD       = FemDataContainer(a);
            mS         = femD.mesh;
            bS         = mS.createBoundaryMesh();
            %             % Generate coordinates
            %             x1 = linspace(0,1,2);
            %             x2 = linspace(0,0.5,2);
            %             % Create the grid
            %             [xv,yv] = meshgrid(x1,x2);
            %             % Triangulate the mesh to obtain coordinates and connectivities
            %             [F,coord] = mesh2tri(xv,yv,zeros(size(xv)),'x');
            %
            %             s.coord    = coord(:,1:2);
            %             s.connec   = F;
            %             mS         = Mesh(s);
            %             bS         = mS.createBoundaryMesh();

            obj.meshReference = mS;
            obj.interfaceMeshReference = bS;
            obj.ninterfaces = length(bS);
        end

        function createDomainMaterial(obj,j,i)
            %             ngaus = 1;
            %             m = obj.meshSubDomain{i};
            obj.material{j,i} = obj.createMaterial(j,i);
        end

        function BC = createRawBoundaryConditionsDirichlet(obj,j,i)
            library = obj.bcDirichletStepLibrary();
            dirichletBc{1,1} = library.dirichletBottomLeft;
            dirichletBc{1,2} = library.dirichletBottomMiddle;
            dirichletBc{1,3} = library.dirichletBottomRight;

            dirichletBc{2,1} = library.dirichletTopLeft;
            dirichletBc{2,2} = library.dirichletTopMiddle;
            dirichletBc{2,3} = library.dirichletTopRight;

            newmanBc{1,1}=[];
            newmanBc{1,2}=[];
            newmanBc{1,3}.boundaryId=2;
            newmanBc{1,3}.dof{1}=[2];
            newmanBc{1,3}.value{1}=[-1];

            newmanBc{2,1}=[];
            newmanBc{2,2}=[];
            newmanBc{2,3}.boundaryId=2;
            newmanBc{2,3}.dof{1}=[2];
            newmanBc{2,3}.value{1}=[-1];

            bM = obj.meshSubDomain{j,i}.createBoundaryMesh();
            [dirichlet,pointload] = obj.createBc(bM,dirichletBc{j,i},newmanBc{j,i});
            BC.dirichlet=dirichlet;
            BC.pointload=pointload;

        end

        function BC = createRawBoundaryConditionsNeumann(obj,j,i)
            library = obj.bcNeumanntStepLibrary();
            dirichletBc{1,1}.boundaryId=[1];
            dirichletBc{1,1}.dof{1}=[1,2];
            dirichletBc{1,1}.value{1}=[0,0];
            dirichletBc{1,2}=[];
            dirichletBc{1,3}=[];

            dirichletBc{2,1}.boundaryId=[1];
            dirichletBc{2,1}.dof{1}=[1,2];
            dirichletBc{2,1}.value{1}=[0,0];
            dirichletBc{2,2}=[];
            dirichletBc{2,3}=[];

            newmanBc{1,1} = library.neumannBottomLeft;
            newmanBc{1,2} = library.neumannBottomMiddle;
            newmanBc{1,3} = library.neumannBottomRight;

            newmanBc{2,1} = library.neumannTopLeft;
            newmanBc{2,2} = library.neumannTopMiddle;
            newmanBc{2,3} = library.neumannTopRight;

           
            bM = obj.meshSubDomain{i}.createBoundaryMesh();
            [dirichlet,pointload] = obj.createBc(bM,dirichletBc{j,i},newmanBc{j,i});
            BC.dirichlet=dirichlet;
            BC.pointload=pointload;

        end

        function [dirichlet,pointload] = createBc(obj,boundaryMesh,dirchletBc,newmanBc)
            if ~isempty(dirchletBc)
                dirichlet = obj.createBoundaryCondition(boundaryMesh,dirchletBc);
            else
                dirichlet=[];
            end
            if ~isempty(newmanBc)
                pointload = obj.createBoundaryCondition(boundaryMesh,newmanBc);
            else
                pointload=[];
            end
        end

        function cond = createBoundaryCondition(obj,bM,condition)
            nbound = length(condition.boundaryId);
            cond = zeros(1,3);
            for ibound=1:nbound
                %                 ncond  = length(condition.dof(nbound,:));
                ncond  = length(condition.dof{ibound});
                nodeId= reshape(unique(bM{condition.boundaryId(ibound)}.globalConnec),[],1);
                nbd   = length(nodeId);
                condition.value{ibound} = condition.value{ibound}/nbd;
                for icond=1:ncond
                    bdcond= [nodeId, repmat(condition.dof{ibound}(icond),[nbd,1]), repmat(condition.value{ibound}(icond),[nbd,1])];
                    cond=[cond;bdcond];
                end
            end
            cond = cond(2:end,:);
        end

        function bc = createBoundaryConditions(obj,mesh,bcV,j,i)
            dim = obj.getFunDims(obj.displacementFun{j,i});
            bcV.ndimf = dim.ndimf;
            bcV.ndofs = dim.ndofs;
            s.mesh  = mesh;
            s.scale = 'MACRO';
            s.bc    = {bcV};
            s.ndofs = dim.ndofs;
            bc = BoundaryConditions(s);
            bc.compute();
        end

        function createSubdomainBoundaryConditions(obj)
            for jdom = 1:obj.nSubdomains(2)
                for idom=1:obj.nSubdomains(1)
                    rawBC = obj.createRawBoundaryConditionsDirichlet(jdom,idom);
                    obj.boundaryConditions.dirichletStep{jdom,idom} = obj.createBoundaryConditions(obj.meshSubDomain{jdom,idom},rawBC,jdom,idom);
                end
            end

            for jdom=1:obj.nSubdomains(2)
                for idom=1:obj.nSubdomains(1)
                    rawBC = obj.createRawBoundaryConditionsNeumann(jdom,idom);
                    obj.boundaryConditions.neumannStep{jdom,idom} = obj.createBoundaryConditions(obj.meshSubDomain{jdom,idom},rawBC,jdom,idom);
                end
            end
        end


        function material = createMaterial(obj,j,i)
            I = ones(obj.meshSubDomain{j,i}.nelem,obj.quad{j,i}.ngaus);
            s.ptype = 'ELASTIC';
            s.pdim  = '2D';
            s.nelem = obj.meshSubDomain{j,i}.nelem;
            s.mesh  = obj.meshSubDomain{j,i};
            s.kappa = .9107*I;
            s.mu    = .3446*I;
            mat = Material.create(s);
            mat.compute(s);
            material = mat;
        end

        function createSubdomainDisplacementFun(obj)
            nx = obj.nSubdomains(1);
            ny = obj.nSubdomains(2);
            for i=1:ny
                for j=1:nx
                    obj.displacementFun{i,j} = P1Function.create(obj.meshSubDomain{i,j}, obj.meshSubDomain{i,j}.ndim);
                end
            end

        end

        function createDdomainDisplacementFun(obj)
            obj.displacementFun= P1Function.create(obj.meshDomain, obj.meshDomain.ndim);
        end

        function dim = getFunDims(obj,disp)
            d.ndimf  = disp.ndimf;
            d.nnodes = size(disp.fValues, 1);
            d.ndofs  = d.nnodes*d.ndimf;
            d.nnodeElem = obj.meshDomain.nnodeElem; % should come from interp..
            d.ndofsElem = d.nnodeElem*d.ndimf;
            dim = d;
        end

        function computeStiffnessMatrix(obj,j,i)
            s.type     = 'ElasticStiffnessMatrix';
            s.mesh     = obj.meshSubDomain{j,i};
            s.fun      = obj.displacementFun{j,i};
            s.material = obj.material{j,i};
            lhs = LHSintegrator.create(s);
            obj.LHS{j,i} = lhs.compute();
        end

        function [Fext,RHS] = computeForces(obj,boundaryConditions,material,mesh,disp,LHS)
            s.type = 'Elastic';
            s.scale    = obj.scale;
            s.dim      = obj.getFunDims(disp);
            s.BC       = boundaryConditions;
            s.mesh     = mesh;
            s.material = material;
            %             s.globalConnec = obj.displacementFun{i}.connec;
            s.globalConnec = mesh.connec;
            RHSint = RHSintegrator.create(s);
            rhs = RHSint.compute();
            Fext = rhs;
            R = RHSint.computeReactions(LHS);
            %             obj.variables.fext = rhs + R;
            RHS = rhs+R;
        end

        function computeSubdomainLHS(obj)
            for jdom = 1:obj.nSubdomains(2)
                for idom=1:obj.nSubdomains(1)
                    obj.quad{jdom,idom} = Quadrature.set(obj.meshSubDomain{jdom,idom}.type);
                    obj.quad{jdom,idom}.computeQuadrature('QUADRATIC');
                    obj.createDomainMaterial(jdom,idom);
                    obj.computeStiffnessMatrix(jdom,idom);
                end
            end
        end

        function computeDomainLHS(obj)
            ndimf  = obj.domainFun.ndimf;
            Gmat   = zeros(obj.meshDomain.nnodes*ndimf);
            for jdom = 1:obj.nSubdomains(2)
                for idom = 1:obj.nSubdomains(1)
                    locMat = obj.LHS{jdom,idom};
                    connec = obj.localGlobalDofConnec{jdom,idom};
                    mat    = obj.local2globalMatrix(locMat,connec);
                    Gmat   = Gmat+mat;
                end
            end
            obj.domLHS = Gmat;
        end



        function hExtG = computeHarmonicExtension(obj,vector)
            hExtG = zeros(size(vector,1),1);
            for jdom= 1:obj.nSubdomains(2)
                for idom = 1:obj.nSubdomains(1)
                    connec  = obj.localGlobalDofConnec{jdom,idom};
                    vec         = obj.global2local(vector,connec);
                    dirDof      = obj.boundaryConditions.neumannStep{jdom,idom}.dirichlet;
                    vec(dirDof) = 0;
                    interfaceDof  = obj.dofInterfaceDomain{jdom,idom};
                    restrictedDof = [interfaceDof; dirDof];
                    interiorDof  = obj.interiorDof{jdom,idom};
                    freeDof = setdiff(interiorDof,restrictedDof);
                    Kii          = obj.LHS{jdom,idom}(freeDof,freeDof);
                    KiL          = obj.LHS{jdom,idom}(freeDof,restrictedDof);
                    hExt   = -Kii\KiL;
                    hExt   = hExt*vec(restrictedDof);
                    hExtD(freeDof) = hExt;
                    hExtD(interfaceDof) = vec(interfaceDof)*0.5;
                    hExtD(dirDof)       = 0;
                    %                     hExtD  = obj.scaleInterfaceValues(hExtD)
                    hExtG   = hExtG + obj.local2global(hExtD,connec);
                end
            end
        end

        function coarseSpace = computeCoarseSpace(obj)
            RBbasisl = obj.computeRigidBodyBasis();
            RBbasisl = obj.scaleInterfaceValues(RBbasisl);
            RBbasisg = obj.createGlobalBasis(RBbasisl);

            %             obj.hExt    = obj.computeHarmonicExtension();
            coarseSpace = obj.computeCoarseSpaceBasis2(RBbasisg);
            %              obj.plotfields(coarseSpace)
            %             coarseSpace = obj.computeCoarseSpaceBasisNoHarmonic(RBbasisg);
        end


        function RBbasis = computeRigidBodyBasis(obj)
            fvalues=[1 0 0;0 1 0; 0 0 1];
            nbasis = size(fvalues,2);
            for jdom = 1:obj.nSubdomains(2)
                for idom = 1:obj.nSubdomains(1)
                    ss.mesh=obj.meshSubDomain{jdom,idom};
                    maxCord = max(ss.mesh.coord);
                    minCord = min(ss.mesh.coord);
                    ss.refPoint= (maxCord + minCord)*0.5;
                    for ibasis = 1:nbasis
                        ss.fvalues = fvalues(ibasis,:);
                        a=RigidBodyFunction(ss);
                        p1FUNC = a.project('P1');
                        bValues = p1FUNC.fValues;
                        bValues = reshape(bValues',1,[])';
                        RBbasis{jdom,idom}(:,ibasis)=bValues;
                    end
                end
            end
        end

        function values = scaleInterfaceValues(obj,val)
            values = val;
            nint = size(obj.interfaceDof,3);
            ndom = size(obj.interfaceDof,2);
            w = [obj.weight,1- obj.weight];
            for iint=1:nint
                for idom = 1:ndom
                    dom = obj.interfaceDom(iint,idom);
                    row = ceil(dom/obj.nSubdomains(1));
                    col = dom-(row-1)*obj.nSubdomains(1);
                    dof = obj.interfaceDof(:,idom,iint);
                    nfields = size(values{row,col},2);
                    for ifield = 1:nfields
                        values{row,col}(dof,ifield) = values{row,col}(dof,ifield)*w(idom);
                    end
                end
            end
        end

        function Gvec = local2global(obj,Lvec,locGlobConnec)
            ndimf  = obj.displacementFun{1}.ndimf;
            Gvec   = zeros(obj.meshDomain.nnodes*ndimf,1);
            Gvec(locGlobConnec(:,1)) = Lvec(locGlobConnec(:,2));
        end

        function Lvec = global2local(obj,Gvec,locGlobConnec)
            ndimf  = obj.displacementFun{1}.ndimf;
            Lvec   = zeros(obj.meshReference.nnodes*ndimf,1);
            Lvec(locGlobConnec(:,2)) = Gvec(locGlobConnec(:,1));
        end

        function globBasis = createGlobalBasis(obj,locBasis)
            nbasis = size(locBasis{1},2);
            for jdom = 1:obj.nSubdomains(2)
                for idom = 1:obj.nSubdomains(1)
                    for ibasis = 1:nbasis
                        locVec  = locBasis{jdom,idom}(:,ibasis);
                        connec  = obj.localGlobalDofConnec{jdom,idom};
                        globVec = obj.local2global(locVec,connec);
                        globBasis{jdom,idom}(:,ibasis) = globVec;
                    end
                end
            end
        end

        function Gmat = local2globalMatrix(obj,Lmat,locGlobConnec)
            ndimf  = obj.domainFun.ndimf;
            Gmat   = zeros(obj.meshDomain.nnodes*ndimf);
            Gmat(locGlobConnec(:,1),locGlobConnec(:,1)) = Lmat(locGlobConnec(:,2),locGlobConnec(:,2));
        end

        function cBasis = computeCoarseSpaceBasis(obj,basis)
            for idom = 1:obj.nSubdomains(1)
                epsilon = obj.hExt{idom};
                nbasis = size(basis{idom},2);
                interfaceDof = obj.dofInterfaceDomain{idom};
                interiorDof = obj.interiorDof{idom};
                for ibasis = 1:nbasis
                    interfaceDisp = basis{idom}(interfaceDof,ibasis);
                    interiorDisp = epsilon*interfaceDisp;
                    cBasis{idom}(interiorDof,ibasis)  = interiorDisp;
                    cBasis{idom}(interfaceDof,ibasis) = interfaceDisp;
                end
            end
        end

        function cBasis = computeCoarseSpaceBasis2(obj,basis)
            nx = size(basis,2);
            ny = size(basis,1);

            for jdom = 1:obj.nSubdomains(2)
                for idom = 1:obj.nSubdomains(1)
                    dvalues = basis{jdom,idom};
                    nbasis = size(dvalues,2);
                    for ibasis = 1:nbasis
                        values = dvalues(:,ibasis);
                        cBasis{jdom,idom}(:,ibasis) = obj.computeHarmonicExtension(values);
                    end
                end
            end

            %             cBasis = obj.scaleInterfaceValues(cBasis);

            %             alldof = 1:1:obj.domainFun.nDofs;
            %             for jdom = 1:obj.nSubdomains(2)
            %                 for idom = 1:obj.nSubdomains(1)
            %                     dofR = obj.localGlobalDofConnec(:,1,idom);
            %                     dofF = setdiff(alldof,dofR);
            %                     Kff = obj.domLHS(dofF,dofF);
            %                     Kfr = obj.domLHS(dofF,dofR);
            %                     nbasis = size(basis{jdom,idom},2);
            %                     for ibasis = 1:nbasis
            %                         disp = basis{idom}(dofR,ibasis);
            %                         F   = Kfr*disp;
            %                         u   = Kff\F;
            %                         cBasis{idom}(dofF,ibasis)  = u;
            %                         cBasis{idom}(dofR,ibasis)  = disp;
            %                     end
            %                 end
            %             end
        end

        function cBasis = computeCoarseSpaceBasisNoHarmonic(obj,basis)
            for idom = 1:obj.nSubdomains(1)
                nbasis = size(basis{idom},2);
                for ibasis = 1:nbasis
                    cBasis{idom}(:,ibasis)  = obj.domLHS*basis{idom}(:,ibasis);
                end
            end
            %             cBasis = obj.domLHS*basis;
        end


        function [interfaceDof,interfaceDom] = computeLocalInterfaceDof(obj)
            intConec = reshape(obj.interfaceConnec',2,obj.interfaceMeshReference{1}.mesh.nnodes,[]);
            intConec = permute(intConec,[2 1 3]);
            nint = size(intConec,3);
            globaldof=0;
            dim = obj.getFunDims(obj.displacementFun{1});
            for iint=1:nint
                ndom = size(intConec,2); %length(intConec(1,:,iint));
                for idom = 1:ndom
                    dofaux=0;
                    nodesI = intConec(:,idom,iint);
                    dom = ceil(intConec(1,idom,iint)/obj.meshReference.nnodes);
                    globaldof = (dom-1)*dim.ndofs;
                    for iunkn=1:dim.ndimf
                        DOF = dim.ndimf*(nodesI - 1) + iunkn;
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

        function dofInterfaceDomain = assingInterfaceDof2Domain(obj)
            nint = size(obj.interfaceDof,3);
            ndom = size(obj.interfaceDof,2);
            dofInterfaceDomain = cell(obj.nSubdomains(2),obj.nSubdomains(1));
            for iint=1:nint
                for idom = 1:ndom
                    dom = obj.interfaceDom(iint,idom);
                    row = ceil(dom/obj.nSubdomains(1));
                    col = dom-(row-1)*obj.nSubdomains(1);
                    dof = obj.interfaceDof(:,idom,iint);
                    if isempty(dofInterfaceDomain{row,col})
                        dofInterfaceDomain{row,col} = dof;
                    else
                        dofInterfaceDomain{row,col} = [dofInterfaceDomain{row,col};dof];
                    end
                end
            end

        end

        function interiorDof = assingInteriorDof(obj)
            for jdom = 1:obj.nSubdomains(2)
                for idom=1:obj.nSubdomains(1)
                    ndof = obj.displacementFun{jdom,idom}.nDofs;
                    dofs = 1:1:ndof;
                    interfaceDof = obj.dofInterfaceDomain{jdom,idom};
                    interiorDof{jdom,idom}  = setdiff(dofs,interfaceDof);
                end
            end
        end

        function  localGlobalDofConnec = createlocalGlobalDofConnec(obj)
            ndimf = obj.displacementFun{1}.ndimf;
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

        function ind = identifyNeumanInterfaceDof(obj)
            for idom=1:obj.nSubdomains(1)
                bN = obj.boundaryConditions.neumannStep{idom};
                for idof = 1: length(obj.interfaceDof(:,idom))
                    ind(idof,idom) = find(bN.neumann == obj.interfaceDof(idof,idom));
                end
            end
        end

        function ind = identifyDirichletInterfaceDof(obj)
            for idom=1:obj.nSubdomains(1)
                bD = obj.boundaryConditions.dirichletStep{idom};
                for idof = 1: length(obj.interfaceDof(:,idom))
                    ind(idof,idom) = find(bD.dirichlet == obj.interfaceDof(idof,idom));
                end
            end
        end


        function NeumannNeumannSolver(obj)
            tol=1e-8;
            e=1;
            iter=1;
            uInt = zeros(size(obj.interfaceDof,1),1);
            while e(iter)>tol

                [uD,RHS] = obj.solveFEM(obj.boundaryConditions.dirichletStep);
                R = obj.computeInterfaceResidual(uD,RHS);
                obj.updateNeumanValues(R);

                [uN,~] = obj.solveFEM(obj.boundaryConditions.neumannStep);
                uInt_new = obj.updateDirichletStepBC(uN);

                iter=iter+1;
                e(iter) = norm(uInt_new-uInt);
                uInt = uInt_new;

                %                 [uD,RHS] = obj.dirichletStep();

                % for idom = 1:obj.nSubdomains(1)
                %      uplot = reshape(obj.displacementFun{idom}.fValues',1,[])';
                %                     uplot = obj.boundaryConditions.neumannStep{idom}.reducedToFullVector(uN{idom});
                %      obj.plotSolution(uplot,obj.meshSubDomain{idom}, idom,iter)
                %  end
                %                 iter=iter+1;

                %                 uN = obj.neumannStep();


                %  for idom = 1:obj.nSubdomains(1)
                %                     uplot = obj.boundaryConditions.neumannStep{idom}.reducedToFullVector(uN{idom});
                %      uplot = uN{idom};
                %      obj.plotSolution(uplot,obj.meshSubDomain{idom}, idom,iter)
                %  end

                %                 iter=iter+1;
                %                 e(iter) = norm(uInt_new-uInt);
                %                 uInt = uInt_new;

                %                iter=iter+1;
                %                 end
            end
        end

        function [uD,RHS] = dirichletStep(obj)
            for idom = 1:obj.nSubdomains(1)
                bD = obj.boundaryConditions.dirichletStep{idom};
                lhs = obj.LHS{idom};
                mesh = obj.meshSubDomain{idom};
                disp = obj.displacementFun{idom};
                mat  = obj.material{idom};
                [Fext,RHS{idom}] = obj.computeForces(bD,mat,mesh,disp,lhs);
                Kred    = bD.fullToReducedMatrix(obj.LHS{idom});
                Fred    = bD.fullToReducedVector(RHS{idom});
                u{idom} = Kred\Fred;
                u{idom} = bD.reducedToFullVector(u{idom});
                uD{idom} = reshape(u{idom},2,[])';
            end
        end

        function u = neumannStep(obj)
            for idom = 1:obj.nSubdomains(1)
                bN = obj.boundaryConditions.neumannStep{idom};
                lhs = obj.LHS{idom};
                mesh = obj.meshSubDomain{idom};
                mat  = obj.material{idom};
                disp = obj.displacementFun{idom};
                [obj.Fext,obj.RHS{idom}] = obj.computeForces(bN,mat,mesh,disp,lhs);
                Kred    = bN.fullToReducedMatrix(obj.LHS{idom});
                Fred    = bN.fullToReducedVector(obj.RHS{idom});
                u{idom} = Kred\Fred;
                u{idom} = bN.reducedToFullVector(u{idom});
            end
        end

        function [u,RHS] = solveFEM(obj,bc)
            for idom = 1:obj.nSubdomains(1)
                %                 bc = obj.boundaryConditions.neumannStep{idom};
                lhs = obj.LHS{idom};
                mesh = obj.meshSubDomain{idom};
                mat  = obj.material{idom};
                disp = obj.displacementFun{idom};
                bc_dom = bc{idom};
                [obj.Fext,RHS{idom}] = obj.computeForces(bc_dom,mat,mesh,disp,lhs);
                Kred    = bc_dom.fullToReducedMatrix(obj.LHS{idom});
                Fred    = bc_dom.fullToReducedVector(RHS{idom});
                u{idom} = Kred\Fred;
                u{idom} = bc_dom.reducedToFullVector(u{idom});
                u{idom} = reshape(u{idom},2,[])';
            end
        end

        function R2 = computeInterfaceResidual(obj,u,RHS)

            w = [obj.weight,1-obj.weight];
            %              w  = [1,1];
            R=zeros(size(obj.interfaceDof,1),1);
            R2=zeros(size(obj.interfaceDof,1),1);
            for idom=1:obj.nSubdomains(1)
                unodal = reshape(u{idom}',1,[])';
                interfaceDOF = obj.interfaceDof(:,idom);
                R2 = R2 + w(idom)*(RHS{idom}(interfaceDOF)-obj.LHS{idom}(interfaceDOF,:)*unodal ...
                    + obj.LHS{idom}(interfaceDOF,interfaceDOF)*unodal(interfaceDOF));
                %                R = R + w(idom)*(obj.Fext{idom}(interfaceDOF)-obj.LHS{idom}(interfaceDOF,:)*unodal);
                %                 aaa = norm(R2-R);
            end
        end

        function updateNeumanValues(obj,R)
            for idom=1:obj.nSubdomains(1)
                ind = obj.interfaceNeumanDof(:,idom);
                obj.boundaryConditions.neumannStep{idom}.neumann_values(ind) = R;
            end
        end


        function uIntNew = updateDirichletStepBC(obj,u)
            uInt = obj.computeInterfaceDisp(u);
            uIntNew = obj.updateDirichletValues(uInt);
        end

        function uInt = computeInterfaceDisp(obj,u)
            w = [obj.weight,1- obj.weight];
            %             w  = [1,1];
            uInt=zeros(size(obj.interfaceDof,1),1);
            for idom=1:obj.nSubdomains(1)
                unodal = reshape(u{idom}',1,[])';
                interfaceDOF = obj.interfaceDof(:,idom);
                uInt = uInt + w(idom)*unodal(interfaceDOF);
            end
        end

        function uIntNew = updateDirichletValues(obj,u)
            for idom=1:obj.nSubdomains(1)
                bD{idom} = obj.boundaryConditions.dirichletStep{idom};
                ind = obj.interfaceDirichletDof(:,idom);
                bD{idom}.dirichlet_values(ind) = bD{idom}.dirichlet_values(ind) + obj.theta*u;

            end
            obj.boundaryConditions.dirichletStep = bD;
            uIntNew = bD{end}.dirichlet_values(obj.interfaceDirichletDof(:,end));
        end

        function plotSolution(obj,x,mesh,domain,iter)
            %             xFull = bc.reducedToFullVector(x);
            s.fValues = reshape(x,2,[])';
            s.mesh = mesh;
            s.fValues(:,end+1) = 0;
            s.ndimf = 3;
            xF = P1Function(s);
            xF.plot();
            %             xF.print(['domain',num2str(domain),'_',num2str(iter)],'Paraview')
            fclose('all');
        end

        function print(obj, filename, software,funNames,fun,mesh)
            if nargin == 2; software = 'GiD'; end
            %             [fun, funNames] = obj.getFunsToPlot();
            a.mesh     = mesh;
            a.filename = filename;
            a.fun      = fun;
            a.funNames = funNames;

            a.type     = software;
            pst = FunctionPrinter.create(a);
            pst.print();
        end

        function plotfields(obj,basis)
            software = 'Paraview';
            funNames = {'translationX', 'translationY', 'rotation'};
            for jdom=1:obj.nSubdomains(2)
                for idom=1:obj.nSubdomains(1)
                    nbasis = size(basis{jdom,idom},2);
                    for ibasis = 1:nbasis
                        values = basis{jdom,idom}(:,ibasis);
                        s.fValues = reshape(values,2,[])';
                        s.mesh = obj.meshDomain;
                        s.fValues(:,end+1) = 0;
                        s.ndimf = 3;
                        fun{ibasis} = P1Function(s);
                    end
                    %                 a.fun=fun;
                    mesh=obj.meshDomain;
                    filename = ['domain',num2str(jdom),num2str(idom)];
                    obj.print(filename,software,funNames,fun,mesh)
                end
            end
        end

        function bc = bcDirichletStepLibrary(obj)
            bc.dirichletBottomLeft.boundaryId=[1,2,4];
            bc.dirichletBottomLeft.dof{1}=[1,2];
            bc.dirichletBottomLeft.dof{2}=[1,2];
            bc.dirichletBottomLeft.dof{3}=[1,2];
            bc.dirichletBottomLeft.value{1}=[0,0];
            bc.dirichletBottomLeft.value{2}=[0,0];
            bc.dirichletBottomLeft.value{3}=[0,0];
            %             bc.newmanNo=[];

            bc.dirichletBottomMiddle.boundaryId=[1,2,4];
            bc.dirichletBottomMiddle.dof{1}=[1,2];
            bc.dirichletBottomMiddle.dof{2}=[1,2];
            bc.dirichletBottomMiddle.dof{3}=[1,2];
            bc.dirichletBottomMiddle.value{1}=[0,0];
            bc.dirichletBottomMiddle.value{2}=[0,0];
            bc.dirichletBottomMiddle.value{3}=[0,0];

            bc.dirichletBottomRight.boundaryId=[1,4];
            bc.dirichletBottomRight.dof{1}=[1,2];
            bc.dirichletBottomRight.dof{2}=[1,2];
            bc.dirichletBottomRight.value{1}=[0,0];
            bc.dirichletBottomRight.value{2}=[0,0];

            bc.dirichletTopLeft.boundaryId=[1,2,3];
            bc.dirichletTopLeft.dof{1}=[1,2];
            bc.dirichletTopLeft.dof{2}=[1,2];
            bc.dirichletTopLeft.dof{3}=[1,2];
            bc.dirichletTopLeft.value{1}=[0,0];
            bc.dirichletTopLeft.value{2}=[0,0];
            bc.dirichletTopLeft.value{3}=[0,0];
            %             bc.newmanNo=[];

            bc.dirichletTopMiddle.boundaryId=[1,2,3];
            bc.dirichletTopMiddle.dof{1}=[1,2];
            bc.dirichletTopMiddle.dof{2}=[1,2];
            bc.dirichletTopMiddle.dof{3}=[1,2];
            bc.dirichletTopMiddle.value{1}=[0,0];
            bc.dirichletTopMiddle.value{2}=[0,0];
            bc.dirichletTopMiddle.value{3}=[0,0];

            bc.dirichletTopRight.boundaryId=[1,3];
            bc.dirichletTopRight.dof{1}=[1,2];
            bc.dirichletTopRight.dof{2}=[1,2];
            %             dirichletBc{3}.dof{2}=[1,2];
            bc.dirichletTopRight.value{1}=[0,0];
            bc.dirichletTopRight.value{2}=[0,0];
        end

        function bc = bcNeumanntStepLibrary(obj)
            bc.neumannBottomLeft.boundaryId=[2,4];
            bc.neumannBottomLeft.dof{1}=[1,2];
            bc.neumannBottomLeft.dof{2}=[1,2];
            bc.neumannBottomLeft.value{1}=[0,0];
            bc.neumannBottomLeft.value{2}=[0,0];
            %             bc.newmanNo=[];

            bc.neumannBottomMiddle.boundaryId=[1,2,4];
            bc.neumannBottomMiddle.dof{1}=[1,2];
            bc.neumannBottomMiddle.dof{2}=[1,2];
            bc.neumannBottomMiddle.dof{3}=[1,2];
            bc.neumannBottomMiddle.value{1}=[0,0];
            bc.neumannBottomMiddle.value{2}=[0,0];
            bc.neumannBottomMiddle.value{3}=[0,0];

            bc.neumannBottomRight.boundaryId=[1,4];
            bc.neumannBottomRight.dof{1}=[1,2];
            bc.neumannBottomRight.dof{2}=[1,2];
            bc.neumannBottomRight.value{1}=[0,0];
            bc.neumannBottomRight.value{2}=[0,0];

            bc.neumannTopLeft.boundaryId=[2,3];
            bc.neumannTopLeft.dof{1}=[1,2];
            bc.neumannTopLeft.dof{2}=[1,2];
            bc.neumannTopLeft.value{1}=[0,0];
            bc.neumannTopLeft.value{2}=[0,0];
            %             bc.newmanNo=[];

            bc.neumannTopMiddle.boundaryId=[1,2,3];
            bc.neumannTopMiddle.dof{1}=[1,2];
            bc.neumannTopMiddle.dof{2}=[1,2];
            bc.neumannTopMiddle.dof{3}=[1,2];
            bc.neumannTopMiddle.value{1}=[0,0];
            bc.neumannTopMiddle.value{2}=[0,0];
            bc.neumannTopMiddle.value{3}=[0,0];

            bc.neumannTopRight.boundaryId=[1,3];
            bc.neumannTopRight.dof{1}=[1,2];
            bc.neumannTopRight.dof{2}=[1,2];
            %             dirichletBc{3}.dof{2}=[1,2];
            bc.neumannTopRight.value{1}=[0,0];
            bc.neumannTopRight.value{2}=[0,0];
        end
    end

    
end
