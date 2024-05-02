classdef BalancingNeumannNeumann_ndom2 < handle

    properties (Access = public)

    end

    properties (Access = private)

    end

    properties (Access = private)
        ndimf
        functionType
        solverCase
        meshDomain
        boundaryConditions
        bcApplier
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
        domainLHS
        LHScoarse

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

        function obj = BalancingNeumannNeumann_ndom2()
            close all
            obj.init();
            obj.createReferenceMesh();

            s.nsubdomains = obj.nSubdomains; %nx ny
            s.meshReference = obj.meshReference;
            m = MeshCreatorFromRVE(s);
            [obj.meshDomain,obj.meshSubDomain,obj.interfaceConnec,~,obj.locGlobConnec] = m.create();

            obj.domainFun = LagrangianFunction.create(obj.meshDomain, obj.ndimf,obj.functionType);
            obj.createSubdomainDisplacementFun();
            [obj.interfaceDof,obj.interfaceDom] = obj.computeLocalInterfaceDof();
            obj.dofInterfaceDomain = obj.assingInterfaceDof2Domain();
            obj.interiorDof        = obj.assingInteriorDof();
            obj.computeSubdomainLHS();
            obj.createSubdomainBoundaryConditions();
            obj.createSubdomainBoundaryConditionsApplier();

            obj.localGlobalDofConnec = obj.createlocalGlobalDofConnec();

            obj.computeDomainLHS();

            [obj.coarseSpace,lambda] = obj.computeCoarseSpace();
%             obj.coarseSpace = obj.coarseSpace(1,1:3);
            obj.LHScoarse = obj.computeCoarseLHS();
            
%             obj.plotfields(lambda)
%             obj.plotfields(obj.coarseSpace)
            

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
            obj.ndimf    = 2;
            obj.nSubdomains = [4 1]; %nx ny
            obj.scale    = 'MACRO';
            obj.weight   = 0.5;
            obj.theta    = 0.2;
            obj.functionType = 'P1';
            obj.solverCase='REDUCED';
        end

        function createReferenceMesh(obj)
%             filename   = 'lattice_ex1';
%             a.fileName = filename;
%             femD       = FemDataContainer(a);
%             mS         = femD.mesh;
%             bS         = mS.createBoundaryMesh();
                        % Generate coordinates
                        x1 = linspace(0,1,5);
                        x2 = linspace(0,0.5,5);
                        % Create the grid
                        [xv,yv] = meshgrid(x1,x2);
                        % Triangulate the mesh to obtain coordinates and connectivities
                        [F,coord] = mesh2tri(xv,yv,zeros(size(xv)),'x');

%                         remove corners
%                         coord(coord(:,1) == max(coord(:,1)) & coord(:,2) == max(coord(:,2)),:)=[];
%                         coord(coord(:,1) == max(coord(:,1)) & coord(:,2) == min(coord(:,2)),:)=[];
%                         coord(coord(:,1) == min(coord(:,1)) & coord(:,2) == max(coord(:,2)),:)=[];
%                         coord(coord(:,1) == min(coord(:,1)) & coord(:,2) == min(coord(:,2)),:)=[];
%                         F = delaunay(coord(:,1:2));
%             
                        s.coord    = coord(:,1:2);
                        s.connec   = F;
                        mS         = Mesh.create(s);
                        bS         = mS.createBoundaryMesh();

            obj.meshReference = mS;
            obj.interfaceMeshReference = bS;
            obj.ninterfaces = length(bS);
        end

        function createDomainMaterial(obj,j,i)
            %             ngaus = 1;
            %             m = obj.meshSubDomain{i};
            obj.material{j,i} = obj.createMaterial(j,i);
        end

        function bc = createRawBoundaryConditionsDirichlet(obj,j,i)
            library = obj.bcDirichletStepLibrary();
            mesh = obj.meshSubDomain{j,i};
            if j==1                   
                if i==1                      
                    dirichletBc = library.dirichletBottomLeft;
                    dirichlet   = DirichletCondition(mesh,dirichletBc);
                    pointload   = [];

                elseif i==obj.nSubdomains(1)
                    dirichletBc = library.dirichletBottomRight;
                    dirichlet   = DirichletCondition(mesh,dirichletBc);

                    isRight      = @(coor) (abs(coor(:,1) - max(coor(:,1)))   < 1e-12);
                    PL.domain    = @(coor) isRight(coor);
                    PL.direction = 2;
                    PL.value     = -1;
                    pointload    = PointLoad(mesh,PL);

                    % need this because force applied in the face not in a point
                    pointload.values=pointload.values/size(pointload.dofs,1);
                    fvalues = zeros(mesh.nnodes*obj.ndimf,1);
                    fvalues(pointload.dofs) = pointload.values;
                    fvalues = reshape(fvalues,obj.ndimf,[])';
                    pointload.fun.fValues = fvalues;

%                     neumannBc{j,i}.boundaryId = 2;
%                     neumannBc{j,i}.dof{1}     = [1,2];
%                     neumannBc{j,i}.value{1}   = [0,-0.1/obj.nSubdomains(1)];
                else                         
                    dirichletBc = library.dirichletBottomMiddle;
                    dirichlet   = DirichletCondition(mesh,dirichletBc);
                    pointload   = [];

                end
            elseif j==obj.nSubdomains(2) 
                if i==1 
                    dirichletBc = library.dirichletTopLeft;
                    dirichlet   = DirichletCondition(mesh,dirichletBc);
                    pointload   = [];
                elseif i==obj.nSubdomains(1)  
                    dirichletBc = library.dirichletTopRight;
                    dirichlet   = DirichletCondition(mesh,dirichletBc);
                    
                    isRight      = @(coor) (abs(coor(:,1) - max(coor(:,1)))   < 1e-12);
                    PL.domain    = @(coor) isRight(coor);
                    PL.direction = 2;
                    PL.value     = -1;
                    pointload    = PointLoad(mesh,PL);

                    % need this because force applied in the face not in a point
                    pointload.values=pointload.values/size(pointload.dofs,1);
                    fvalues = zeros(mesh.nnodes*obj.ndimf,1);
                    fvalues(pointload.dofs) = pointload.values;
                    fvalues = reshape(fvalues,obj.ndimf,[])';
                    pointload.fun.fValues = fvalues;
                else      
                    dirichletBc = library.dirichletTopMiddle;
                    dirichlet   = DirichletCondition(mesh,dirichletBc);
                    pointload   = [];
                end
            elseif i==1
                dirichletBc = library.dirichletLeftMiddle;
                dirichlet   = DirichletCondition(mesh,dirichletBc);
                pointload   = [];

            elseif i==obj.nSubdomains(1)
                dirichletBc = library.dirichletRightMiddle;
                dirichlet   = DirichletCondition(mesh,dirichletBc);
                isRight      = @(coor) (abs(coor(:,1) - max(coor(:,1)))   < 1e-12);
                PL.domain    = @(coor) isRight(coor);
                PL.direction = 2;
                PL.value     = -1;
                pointload    = PointLoad(mesh,PL);

                % need this because force applied in the face not in a point
                pointload.values=pointload.values/size(pointload.dofs,1);
                fvalues = zeros(mesh.nnodes*obj.ndimf,1);
                fvalues(pointload.dofs) = pointload.values;
                fvalues = reshape(fvalues,obj.ndimf,[])';
                pointload.fun.fValues = fvalues;
            else
                dirichletBc = library.dirichletMiddleMiddle;
                dirichlet   = DirichletCondition(mesh,dirichletBc);
                pointload   = [];
            end

            s.pointloadFun = pointload;
            s.dirichletFun = dirichlet;
            s.periodicFun  = [];
            s.mesh         = mesh;
            bc             = BoundaryConditions(s);
%             bM = obj.meshSubDomain{j,i}.createBoundaryMesh();
%             [dirichlet,pointload] = obj.createBc(bM,dirichletBc{j,i},neumannBc{j,i});
%             BC.dirichlet=dirichlet;
%             BC.pointload=pointload;

        end

        function bc = createRawBoundaryConditionsNeumann(obj,j,i)
            library = obj.bcNeumannStepLibrary();
            mesh = obj.meshSubDomain{j,i};
            if j==1
                if i==1
                    newmanBc    = library.neumannBottomLeft;
                    pointload   = PointLoad(mesh,newmanBc);

                    isLeft                = @(coor) (abs(coor(:,1) - min(coor(:,1)))   < 1e-12);
                    dirichletBc.domain    = @(coor) isLeft(coor);
                    dirichletBc.direction = [1,2];
                    dirichletBc.value     = 0;

                    dirichlet = DirichletCondition(mesh,dirichletBc);
%                     dirichletBc{j,i}.boundaryId = [1];
%                     dirichletBc{j,i}.dof{1}     = [1,2];
%                     dirichletBc{j,i}.value{1}   = [0,0];
                elseif i==obj.nSubdomains(1)
                    newmanBc  = library.neumannBottomRight;
                    pointload = PointLoad(mesh,newmanBc);
                    dirichlet = [];
                else
                    newmanBc  = library.neumannBottomMiddle;
                    pointload = PointLoad(mesh,newmanBc);
                    dirichlet = [];
                end
            elseif j==obj.nSubdomains(2)
                if i==1
                    newmanBc  = library.neumannTopLeft;
                    pointload = PointLoad(mesh,newmanBc);

                    isLeft                = @(coor) (abs(coor(:,1) - min(coor(:,1)))   < 1e-12);
                    dirichletBc.domain    = @(coor) isLeft(coor);
                    dirichletBc.direction = [1,2];
                    dirichletBc.value     = 0;

                    dirichlet = DirichletCondition(mesh,dirichletBc);
                elseif i==obj.nSubdomains(1)
                    newmanBc  = library.neumannTopRight;
                    pointload = PointLoad(mesh,newmanBc);
                    dirichlet = [];
                else
                    newmanBc  = library.neumannTopMiddle;
                    pointload = PointLoad(mesh,newmanBc);
                    dirichlet = [];
                end
            elseif i==1
                newmanBc  = library.neumannLeftMiddle;
                pointload = PointLoad(mesh,newmanBc);

                isLeft                = @(coor) (abs(coor(:,1) - min(coor(:,1)))   < 1e-12);
                dirichletBc.domain    = @(coor) isLeft(coor);
                dirichletBc.direction = [1,2];
                dirichletBc.value     = 0;
                dirichlet = DirichletCondition(mesh,dirichletBc);
            elseif i==obj.nSubdomains(1)
                newmanBc  = library.neumannRightMiddle;
                pointload = PointLoad(mesh,newmanBc);
                dirichlet = [];
            else
                newmanBc  = library.neumannMiddleMiddle;
                pointload = PointLoad(mesh,newmanBc);
                dirichlet = [];
            end

            s.pointloadFun = pointload;
            s.dirichletFun = dirichlet;
            s.periodicFun  = [];
            s.mesh         = mesh;
            bc             = BoundaryConditions(s);

%             bM = obj.meshSubDomain{i}.createBoundaryMesh();
%             [dirichlet,pointload] = obj.createBc(bM,dirichletBc{j,i},newmanBc{j,i});
%             BC.dirichlet=dirichlet;
%             BC.pointload=pointload;

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
                     obj.boundaryConditions.dirichletStep{jdom,idom} = rawBC;
%                     obj.boundaryConditions.dirichletStep{jdom,idom} = obj.createBoundaryConditions(obj.meshSubDomain{jdom,idom},rawBC,jdom,idom);
                end
            end

            for jdom=1:obj.nSubdomains(2)
                for idom=1:obj.nSubdomains(1)
                    rawBC = obj.createRawBoundaryConditionsNeumann(jdom,idom);
                    obj.boundaryConditions.neumannStep{jdom,idom} = rawBC;
%                     obj.boundaryConditions.neumannStep{jdom,idom} = obj.createBoundaryConditions(obj.meshSubDomain{jdom,idom},rawBC,jdom,idom);
                end
            end
        end

        function createSubdomainBoundaryConditionsApplier(obj)
            for jdom = 1:obj.nSubdomains(2)
                for idom=1:obj.nSubdomains(1)
                    s.boundaryConditions   = obj.boundaryConditions.dirichletStep{jdom,idom};
                    s.mesh                 = obj.meshSubDomain{jdom,idom};
                    obj.bcApplier.dirichletStep{jdom,idom} = BCApplier(s);
                end
            end

            for jdom = 1:obj.nSubdomains(2)
                for idom=1:obj.nSubdomains(1)
                    s.boundaryConditions   = obj.boundaryConditions.neumannStep{jdom,idom};
                    s.mesh                 = obj.meshSubDomain{jdom,idom};
                    obj.bcApplier.neumannStep{jdom,idom} = BCApplier(s);
                end
            end

        end

        function [young,poisson] = computeElasticProperties(obj,mesh)
            E1  = 1;
            nu1 = 1/3;
            E   = AnalyticalFunction.create(@(x) E1*ones(size(squeeze(x(1,:,:)))),1,mesh);
            nu  = AnalyticalFunction.create(@(x) nu1*ones(size(squeeze(x(1,:,:)))),1,mesh);
            young   = E;
            poisson = nu;
        end

        function material = createMaterial(obj,j,i)
            mesh = obj.meshSubDomain{j,i};
            [young,poisson] = obj.computeElasticProperties(mesh);
            s.type    = 'ISOTROPIC';
            s.ptype   = 'ELASTIC';
            s.ndim    = mesh.ndim;
            s.young   = young;
            s.poisson = poisson;
            tensor    = Material.create(s);
            material  = tensor;
        end

%         function material = createMaterial(obj,j,i)
%             I = ones(obj.meshSubDomain{j,i}.nelem,obj.quad{j,i}.ngaus);
%             s.ptype = 'ELASTIC';
%             s.pdim  = '2D';
%             s.nelem = obj.meshSubDomain{j,i}.nelem;
%             s.mesh  = obj.meshSubDomain{j,i};
%             s.kappa = .9107*I;
%             s.mu    = .3446*I;
%             mat = Material.create(s);
%             mat.compute(s);
%             material = mat;
%         end

        function createSubdomainDisplacementFun(obj)
            nx = obj.nSubdomains(1);
            ny = obj.nSubdomains(2);
            for i=1:ny
                for j=1:nx
                    obj.displacementFun{i,j} = LagrangianFunction.create(obj.meshSubDomain{i,j},obj.ndimf,obj.functionType);
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

%         function computeStiffnessMatrix(obj,j,i)
%             s.type     = 'ElasticStiffnessMatrix';
%             s.mesh     = obj.meshSubDomain{j,i};
%             s.fun      = obj.displacementFun{j,i};
%             s.material = obj.material{j,i};
%             lhs = LHSintegrator.create(s);
%             obj.LHS{j,i} = lhs.compute();
%         end

        function LHS = computeStiffnessMatrix(obj,j,i)
            s.type     = 'ElasticStiffnessMatrix';
            s.mesh     = obj.meshSubDomain{j,i};
            s.test     = LagrangianFunction.create(s.mesh,obj.ndimf, 'P1');
            s.trial    = obj.displacementFun{j,i};
            s.material = obj.material{j,i};
            s.quadratureOrder = 2;
            lhs = LHSintegrator.create(s);
            LHS = lhs.compute();
        end

%         function [Fext,RHS] = computeForces(obj,boundaryConditions,material,mesh,disp,LHS)
%             s.type = 'Elastic';
%             s.scale    = obj.scale;
%             s.dim      = obj.getFunDims(disp);
%             s.BC       = boundaryConditions;
%             s.mesh     = mesh;
%             s.material = material;
%             %             s.globalConnec = obj.displacementFun{i}.connec;
%             s.globalConnec = mesh.connec;
%             RHSint = RHSintegrator.create(s);
%             rhs = RHSint.compute();
%             Fext = rhs;
%             R = RHSint.computeReactions(LHS);
%             %             obj.variables.fext = rhs + R;
%             RHS = rhs+R;
%         end

         function forces = computeForces(obj,mesh,dispFun,boundaryConditions,material,stiffness)
            s.type      = 'Elastic';
            s.scale     = 'MACRO';
%             s.dim       = obj.getFunDims();
            s.dim.ndofs = dispFun.nDofs;
            s.BC       = boundaryConditions;
            s.mesh     = mesh;
            s.material = material;
%             s.globalConnec = mesh.connec;
            RHSint = RHSintegrator.create(s);
            rhs = RHSint.compute();
            % Perhaps move it inside RHSint?
            if strcmp(obj.solverCase,'REDUCED')
                R = RHSint.computeReactions(stiffness);
                forces = rhs+R;
            else
                forces = rhs;
            end
        end


        function computeSubdomainLHS(obj)
            for jdom = 1:obj.nSubdomains(2)
                for idom=1:obj.nSubdomains(1)
%                     obj.quad{jdom,idom} = Quadrature.set(obj.meshSubDomain{jdom,idom}.type);
%                     obj.quad{jdom,idom}.computeQuadrature('QUADRATIC');
                    obj.createDomainMaterial(jdom,idom);
                    obj.LHS{jdom,idom} = obj.computeStiffnessMatrix(jdom,idom);
                end
            end
        end

        function computeDomainLHS(obj)
            Gmat   = zeros(obj.meshDomain.nnodes*obj.ndimf);
            for jdom = 1:obj.nSubdomains(2)
                for idom = 1:obj.nSubdomains(1)
                    locMat = obj.LHS{jdom,idom};
                    connec = obj.localGlobalDofConnec{jdom,idom};
                    mat    = obj.local2globalMatrix(locMat,connec);
                    Gmat   = Gmat+mat;
                end
            end
            obj.domainLHS = Gmat;
        end



        function hExtG = computeHarmonicExtension(obj,vector)
            hExtG = zeros(size(vector,1),1);
            for jdom= 1:obj.nSubdomains(2)
                for idom = 1:obj.nSubdomains(1)
                    connec  = obj.localGlobalDofConnec{jdom,idom};
                    vec         = obj.global2local(vector,connec);
                    dirDof      = obj.boundaryConditions.neumannStep{jdom,idom}.dirichlet_dofs;
                    neuDof      = [obj.boundaryConditions.dirichletStep{jdom,idom}.pointload_dofs];
%                     neuDof      = [];
                    vec(dirDof) = 0;
                    interfaceDof  = obj.dofInterfaceDomain{jdom,idom};
                    restrictedDof = [interfaceDof; dirDof; neuDof];
                    interiorDof  = obj.interiorDof{jdom,idom};
                    freeDof = setdiff(interiorDof,restrictedDof);
                    Kii          = obj.LHS{jdom,idom}(freeDof,freeDof);
                    KiL          = obj.LHS{jdom,idom}(freeDof,restrictedDof);
                    hExt   = -Kii\KiL;
                    hExt   = hExt*vec(restrictedDof);
                    hExtD(freeDof) = hExt;
                    hExtD(interfaceDof) = vec(interfaceDof)*0.5;
                    hExtD(dirDof)       = 0;
                    hExtD(neuDof)       = vec(neuDof);  
                    %                     hExtD  = obj.scaleInterfaceValues(hExtD)
                    hExtG   = hExtG + obj.local2global(hExtD,connec);
                end
            end
        end

        function [coarseSpace,lambda] = computeCoarseSpace(obj)
            RBbasisl = obj.computeRigidBodyBasis();
%             RBbasisl = obj.scaleInterfaceValues(RBbasisl);
%             RBbasisg = obj.createGlobalBasis(RBbasisl);

            %             obj.hExt    = obj.computeHarmonicExtension();
            [coarseSpace,lambda] = obj.computeCoarseSpaceBasis2(RBbasisl);

            
            %              obj.plotfields(coarseSpace)
            %             coarseSpace = obj.computeCoarseSpaceBasisNoHarmonic(RBbasisg);
        end


        function RBbasis = computeRigidBodyBasis(obj)
            fvalues=[1 0 0;0 1 0; 0 0 1];
            nbasis = size(fvalues,2);
            ind=1;
            for jdom = 1:obj.nSubdomains(2)
                for idom = 1:obj.nSubdomains(1)
                    ss.mesh = obj.meshSubDomain{jdom,idom};
                    connec  = obj.localGlobalDofConnec{jdom,idom};
                    maxCord = max(ss.mesh.coord);
                    minCord = min(ss.mesh.coord);
                    ss.refPoint= (maxCord + minCord)*0.5;
                    if isempty(obj.boundaryConditions.neumannStep{jdom,idom}.dirichlet_dofs)
                        for ibasis = 1:nbasis
                            ss.fvalues = fvalues(ibasis,:);
                            a=RigidBodyFunction(ss);
                            p1FUNC = a.project('P1');
                            bValues = p1FUNC.fValues;
                            bValues = reshape(bValues',1,[])';
                            bValues = obj.local2global(bValues,connec);
                            RBbasis{ind}(:,ibasis)=bValues;
                        end
                        ind=ind+1;
                    end
                end
            end
        end

        function values = scaleInterfaceValues(obj,val)
            values = val;
            nint = size(obj.interfaceDof,3);
            ndom = size(obj.interfaceDof,2);
            w = [obj.weight,1- obj.weight];
            nbasis = size(values,2);
            for ibasis = 1:nbasis
                nfields = size(values{ibasis},2);
                for ifield=1:nfields
                    gVec = values{ibasis}(:,ifield);
%                     sVec = zeros(size(gVec,1),1);
                    for iint=1:nint
                         sVec = zeros(size(gVec,1),1);
                        for idom = 1:ndom
                            dom = obj.interfaceDom(iint,idom);
                            row = ceil(dom/obj.nSubdomains(1));
                            col = dom-(row-1)*obj.nSubdomains(1);
                            dof = obj.interfaceDof(:,idom,iint);
                            connec = obj.localGlobalDofConnec{row,col};
                            lValue = obj.global2local(gVec,connec);
                            lValue(dof) = lValue(dof)*w(idom);
                            gsVec = obj.local2global(lValue,connec);
                            sVec = sVec + gsVec;
                        end
                        gVec = sVec;
                    end
                    values{ibasis}(:,ifield) = sVec;
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

        function [cBasis,lambda] = computeCoarseSpaceBasis2(obj,basis)
            nx = size(basis,2);
            ny = size(basis,1);
            nbasis = size(basis,2);
            for ibasis = 1:nbasis
                nfields = size(basis{ibasis},2);
                values  = basis{ibasis};
                for ifield = 1:nfields
                     ivalues = values(:,ifield)*0.5;
                     cBasis{ibasis}(:,ifield) =  obj.computeHarmonicExtension(ivalues);
                     lambda{ibasis}(:,ifield) =  -obj.domainLHS*cBasis{ibasis}(:,ifield);
                end
            end

%             for jdom = 1:obj.nSubdomains(2)
%                 for idom = 1:obj.nSubdomains(1)
%                     dvalues = basis{jdom,idom};
%                     nbasis = size(dvalues,2);
%                     for ibasis = 1:nbasis
%                         values = dvalues(:,ibasis);
%                         cBasis{jdom,idom}(:,ibasis) = obj.computeHarmonicExtension(values);
%                     end
%                 end
%             end
        end

        function cBasis = computeCoarseSpaceBasisNoHarmonic(obj,basis)
            for idom = 1:obj.nSubdomains(1)
                nbasis = size(basis{idom},2);
                for ibasis = 1:nbasis
                    cBasis{idom}(:,ibasis)  = obj.domainLHS*basis{idom}(:,ibasis);
                end
            end
            %             cBasis = obj.domLHS*basis;
        end

        function Kcoarse = computeCoarseLHS(obj)
            Kfine = obj.domainLHS;
            nbasis = size(obj.coarseSpace,2);
            for ibasis = 1:nbasis
                 basis = obj.coarseSpace{ibasis};
                 Kcoarse{ibasis} = obj.projectMat(Kfine,basis);
            end
        end

        function pmat = projectMat(obj,mat,basis)
            pmat = mat*basis;
            pmat = basis'*pmat;
        end

        function pvec = projectVec(obj,vec,basis)
            pvec = basis'*vec;
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
            nint = size(obj.interfaceDof,3);
            for iint = 1:nint
                ndom = size(obj.interfaceDof(:,:,iint),2);
                for idom = 1:ndom
                    dom = obj.interfaceDom(iint,idom);
                    row = ceil(dom/obj.nSubdomains(1));
                    col = dom-(row-1)*obj.nSubdomains(1);
                    dof = obj.interfaceDof(:,idom,iint);
                    bN = obj.boundaryConditions.neumannStep{row,col};
                    [~,ind(:,idom,iint)] = ismember(dof,bN.pointload_dofs);
                end
            end
        end

        function ind = identifyDirichletInterfaceDof(obj)
            nint = size(obj.interfaceDof,3);
            for iint = 1:nint
                ndom = size(obj.interfaceDof(:,:,iint),2);
                for idom = 1:ndom
                    dom = obj.interfaceDom(iint,idom);
                    row = ceil(dom/obj.nSubdomains(1));
                    col = dom-(row-1)*obj.nSubdomains(1);
                    dof = obj.interfaceDof(:,idom,iint);
                    bN = obj.boundaryConditions.dirichletStep{row,col};
                    [~,ind(:,idom,iint)] = ismember(dof,bN.dirichlet_dofs);
                end
            end
%             for idom=1:obj.nSubdomains(1)
%                 bD = obj.boundaryConditions.dirichletStep{idom};
%                 for idof = 1: length(obj.interfaceDof(:,idom))
%                     ind(idof,idom) = find(bD.dirichlet == obj.interfaceDof(idof,idom));
%                 end
%             end
        end


        function NeumannNeumannSolver(obj)
            tol=1e-8;
            e=1;
            iter=1;
            nint = size(obj.interfaceDof,3);
            uIntOld = zeros(size(obj.interfaceDof,1),nint);
            while e(iter)>tol

                [uD,RHS] = obj.solveFEM(obj.boundaryConditions.dirichletStep,obj.bcApplier.dirichletStep);
                R = obj.computeInterfaceResidual(uD,RHS);
%                 obj.updateNeumanValues(R);
                RG = obj.constructGlobalResidual(R);
%                 obj.plotSolution(RG,obj.meshDomain,11,11,iter,2)
%                 for jdom = 1:obj.nSubdomains(2)
%                     for idom =1:obj.nSubdomains(1)
%                         mesh = obj.meshSubDomain{jdom,idom};
%                         x = uD{jdom,idom};
%                         row = jdom;
%                         col = idom;
%                         obj.plotSolution(x,mesh,row,col,iter)
%                     end
%                 end

%                 [u0,RHS] = obj.solveCoarseProblem(obj.boundaryConditions.neumannStep);
                u0 = obj.solveCoarseProblem2(RG);
%                 obj.plotSolution(u0,obj.meshDomain,10,10,iter)
                Rbal = obj.balanceResidual3(u0,RG);
%                 for jdom = 1:obj.nSubdomains(2)
%                     for idom =1:obj.nSubdomains(1)
%                         mesh = obj.meshSubDomain{jdom,idom};
%                         x = Rbal{jdom,idom};
%                         row = jdom;
%                         col = idom;
%                         obj.plotSolution(x,mesh,row,col,iter,1)
%                     end
%                 end
%                 Rbal = obj.balanceResidual2(u0,RHS);
                obj.updateBalancedNeumanValues(Rbal);
%                 R = obj.computeInterfaceResidual(u0,RHS);
%                 obj.updateNeumanValues(R);

                [uN,~] = obj.solveFEM(obj.boundaryConditions.neumannStep,obj.bcApplier.dirichletStep);
%                 for jdom = 1:obj.nSubdomains(2)
%                     for idom =1:obj.nSubdomains(1)
%                         mesh = obj.meshSubDomain{jdom,idom};
%                         x = uN{jdom,idom};
%                         row = jdom;
%                         col = idom;
%                         obj.plotSolution(x,mesh,row,col,iter,3)
%                     end
%                 end
%                 [uN,~] = obj.solveFineNeumann(obj.boundaryConditions.neumannStep,Rbal);
                uN     = obj.constructNeumannDisplacement(uN,u0); 

%                  for jdom = 1:obj.nSubdomains(2)
%                     for idom =1:obj.nSubdomains(1)
%                         mesh = obj.meshSubDomain{jdom,idom};
%                         x = uN{jdom,idom};
%                         row = jdom;
%                         col = idom;
%                         obj.plotSolution(x,mesh,row,col,iter,4)
%                     end
%                 end
%                 uN = u0 + uN;;
                uInt = obj.computeInterfaceDisp(uN);
                uIntNew = obj.updateDirichletValues(uInt);

                iter=iter+1;
                e(iter) = norm(uIntNew-uIntOld)/norm(uIntOld);
                if e(iter)>e(iter-1)
                    aaaa=1;
                end
                uIntOld = uIntNew;

              
            end
        end

        function [uD,RHS] = dirichletStep(obj)
            for jdom=1:obj.nSubdomains(2)
                for idom = 1:obj.nSubdomains(1)
                    bD = obj.boundaryConditions.dirichletStep{jdom,idom};
                    lhs = obj.LHS{jdom,idom};
                    mesh = obj.meshSubDomain{jdom,idom};
                    disp = obj.displacementFun{jdom,idom};
                    mat  = obj.material{jdom,idom};
                    [Fext,RHS{idom}] = obj.computeForces(bD,mat,mesh,disp,lhs);
                    Kred    = bD.fullToReducedMatrix(obj.LHS{idom});
                    Fred    = bD.fullToReducedVector(RHS{idom});
                    u{idom} = Kred\Fred;
                    u{idom} = bD.reducedToFullVector(u{idom});
                    uD{idom} = reshape(u{idom},2,[])';
                end
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

        function [u,RHS] = solveFEM(obj,bc,bcApplier)
            for jdom=1:obj.nSubdomains(2)
                for idom = 1:obj.nSubdomains(1)
                    %                 bc = obj.boundaryConditions.neumannStep{idom};
                    lhs = obj.LHS{jdom,idom};
                    mesh = obj.meshSubDomain{jdom,idom};
                    mat  = obj.material{jdom,idom};
                    disp = obj.displacementFun{jdom,idom};
                    bc_dom = bc{jdom,idom};
                    bcApp = bcApplier{jdom,idom};
                    RHS{jdom,idom} = obj.computeForces(mesh,disp,bc_dom,mat,lhs);
                    Kred    = bcApp.fullToReducedMatrixDirichlet(lhs);
                    Fred    = bcApp.fullToReducedVectorDirichlet(RHS{jdom,idom});
                    uRed    = pinv(full(Kred))*Fred;
                    uF      = bcApp.reducedToFullVectorDirichlet(uRed);
                    u{jdom,idom} = reshape(uF,2,[])';
                end
            end
        end

        function [u,RHS] = solveFineNeumann(obj,bc,RHS)
            for jdom=1:obj.nSubdomains(2)
                for idom = 1:obj.nSubdomains(1)
                    bc_dom = bc{jdom,idom};
                    Kred    = bc_dom.fullToReducedMatrix(obj.LHS{jdom,idom});
                    Fred    = bc_dom.fullToReducedVector(RHS{jdom,idom});
                    uRed = pinv(full(Kred))*Fred;
%                     uRed = Kred\Fred;
                    uF = bc_dom.reducedToFullVector(uRed);
                    u{jdom,idom} = reshape(uF,2,[])';
                end
            end
        end

        function [uG,RHS] = solveCoarseProblem(obj,bc)
            ndimf  = obj.domainFun.ndimf;
            uG   = zeros(obj.meshDomain.nnodes*ndimf,1);
            nbasis = size(obj.coarseSpace,2);
            for ibasis = 1:nbasis
                basis = obj.coarseSpace{ibasis};
                for jdom=1:obj.nSubdomains(2)
                    for idom = 1:obj.nSubdomains(1)
                        %                 bc = obj.boundaryConditions.neumannStep{idom};
                        lhs = obj.LHS{jdom,idom};
                        mesh = obj.meshSubDomain{jdom,idom};
                        mat  = obj.material{jdom,idom};
                        disp = obj.displacementFun{jdom,idom};
                        bc_dom = bc{jdom,idom};
                        RHS{jdom,idom} = obj.computeForces(mesh,disp,bc_dom,mat,lhs);
%                         basis = obj.coarseSpace{jdom,idom};
                        connec = obj.localGlobalDofConnec{jdom,idom};
                        F = obj.local2global(RHS{jdom,idom},connec);
                        F = obj.projectVec(F,basis)*0.5;
                        uRed    = obj.LHScoarse{ibasis}\F;
                        uGb      = obj.projectVec(uRed,basis');
                        %                     uL      = obj.global2local(uG,connec);
                        
                            uG = uG + uGb;
                        
                        %                     uG = uG + uGb;
                        %                     u{jdom,idom} = reshape(uG,2,[])';
                        %                     u{jdom,idom} = uL;
                    end
                end
%                 uG = reshape(uG,2,[])';
            end
            uG = reshape(uG,2,[])';
        end

        function uG = solveCoarseProblem2(obj,R)
            ndimf  = obj.domainFun.ndimf;
            uG   = zeros(obj.meshDomain.nnodes*ndimf,1);
            nbasis = size(obj.coarseSpace,2);
            for ibasis = 1:nbasis
                basis = obj.coarseSpace{ibasis};
                F     = obj.projectVec(R,basis);
                uRed  = obj.LHScoarse{ibasis}\F;
                uGb   = obj.projectVec(uRed,basis');
                uG = uG + uGb;
            end
            uG = reshape(uG,2,[])';
        end

        function Rint = computeInterfaceResidual(obj,u,RHS)

            w = [obj.weight,1-obj.weight];
            %              w  = [1,1];
            R=zeros(size(obj.interfaceDof,1),1);
%             R2=zeros(size(obj.interfaceDof,1),1);
            nint = size(obj.interfaceDof,3);
            for iint = 1:nint
                ndom = size(obj.interfaceDof(:,:,iint),2);
                R=zeros(size(obj.interfaceDof,1),1);
                for idom = 1:ndom
                    dom = obj.interfaceDom(iint,idom);
                    row = ceil(dom/obj.nSubdomains(1));
                    col = dom-(row-1)*obj.nSubdomains(1);
                    unodal = reshape(u{row,col}',1,[])';
                    dof = obj.interfaceDof(:,idom,iint);
%                     R = R + w(idom)*(RHS{row,col}(dof)-obj.LHS{row,col}(dof,:)*unodal ...
%                     + obj.LHS{row,col}(dof,dof)*unodal(dof));
                    R = R + 1*(RHS{row,col}(dof)-obj.LHS{row,col}(dof,:)*unodal ...
                    + obj.LHS{row,col}(dof,dof)*unodal(dof));
                end
                Rint(:,iint)=R;
            end
%             for idom=1:obj.nSubdomains(1)
%                 unodal = reshape(u{idom}',1,[])';
%                 interfaceDOF = obj.interfaceDof(:,idom);
%                 R2 = R2 + w(idom)*(RHS{idom}(interfaceDOF)-obj.LHS{idom}(interfaceDOF,:)*unodal ...
%                     + obj.LHS{idom}(interfaceDOF,interfaceDOF)*unodal(interfaceDOF));
%                 %                R = R + w(idom)*(obj.Fext{idom}(interfaceDOF)-obj.LHS{idom}(interfaceDOF,:)*unodal);
%                 %                 aaa = norm(R2-R);
%             end
        end

        function RG = constructGlobalResidual(obj,interfaceR)
            nint = size(obj.interfaceDof,3);
            ndimf  = obj.domainFun.ndimf;
            RG   = zeros(obj.meshDomain.nnodes*ndimf,1);
            for iint = 1:nint
                dom = obj.interfaceDom(iint,1);
                row = ceil(dom/obj.nSubdomains(1));
                col = dom-(row-1)*obj.nSubdomains(1);
                dof = obj.interfaceDof(:,1,iint);
                dim = obj.getFunDims(obj.displacementFun{row,col});
                ndof = dim.ndofs;
                RL =  zeros(ndof,1);
                RL(dof) = interfaceR(:,iint);
                connec = obj.localGlobalDofConnec{row,col};
                RGi = obj.local2global(RL,connec);
                RG = RG + RGi;
            end
        end

        function updateNeumanValues(obj,R)
            nint = size(obj.interfaceDof,3);
            for iint = 1:nint
                ndom = size(obj.interfaceDof(:,:,iint),2);
                for idom = 1:ndom
                    dom = obj.interfaceDom(iint,idom);
                    row = ceil(dom/obj.nSubdomains(1));
                    col = dom-(row-1)*obj.nSubdomains(1);
                    ind = obj.interfaceNeumanDof(:,idom,iint);
                    obj.boundaryConditions.neumannStep{row,col}.neumann_values(ind) = R(:,iint);
                end
            end
        end

%         function updateBalancedNeumanValues2(obj,R)
%             nint = size(obj.interfaceDof,3);
%             for iint = 1:nint
%                 ndom = size(obj.interfaceDof(:,:,iint),2);
%                 for idom = 1:ndom
%                     dom = obj.interfaceDom(iint,idom);
%                     row = ceil(dom/obj.nSubdomains(1));
%                     col = dom-(row-1)*obj.nSubdomains(1);
%                     dof = obj.interfaceDof(:,idom,iint);
%                     ind = obj.interfaceNeumanDof(:,idom,iint);
%                     obj.boundaryConditions.neumannStep{row,col}.neumann_values(ind) = R{row,col}(dof);
%                 end
%             end
%         end

        function updateBalancedNeumanValues(obj,R)
            nint = size(obj.interfaceDof,3);
            for iint = 1:nint
                ndom = size(obj.interfaceDof(:,:,iint),2);
                for idom = 1:ndom
                    dom = obj.interfaceDom(iint,idom);
                    row = ceil(dom/obj.nSubdomains(1));
                    col = dom-(row-1)*obj.nSubdomains(1);
%                     ind = obj.interfaceNeumanDof(:,idom,iint);
                    dof = obj.interfaceDof(:,idom,iint);
                    bc  = obj.boundaryConditions.neumannStep{row,col};
                    values = reshape(bc.pointloadFun.fValues',1,[])';
                    values(dof) = R{row,col}(dof);
                    values = reshape(values,obj.ndimf,[])';
                    obj.boundaryConditions.neumannStep{row,col}.pointloadFun.fValues = values;
                end
            end
        end

        function u = constructNeumannDisplacement(obj,uFine,uCoarse)          
            uC = reshape(uCoarse',1,[])';
            for jdom = 1:obj.nSubdomains(2)
                for idom = 1:obj.nSubdomains(1)
                    connec = obj.localGlobalDofConnec{jdom,idom};
                    uF = reshape(uFine{jdom,idom}',1,[])';
                    uCL = obj.global2local(uC,connec);
                    uL  = uF + uCL;
                    u{jdom,idom} = reshape(uL,2,[])';
                end
            end
        end

        function Rbal = balanceResidual3(obj,u,RG)
            unodal = reshape(u',1,[])';
            rCoarse = obj.domainLHS*unodal;
            RGbal = RG-rCoarse;
            for jdom=1:obj.nSubdomains(2)
                for idom =1:obj.nSubdomains(1)
                    connec = obj.localGlobalDofConnec{jdom,idom};
                    Rbal{jdom,idom} = obj.global2local(RGbal,connec);
                end
            end          
        end
        
        
        function Rbal = balanceResidual2(obj,u,RHS)
            unodal = reshape(u',1,[])';
            rCoarse = obj.domainLHS*unodal;
            for jdom=1:obj.nSubdomains(2)
                for idom =1:obj.nSubdomains(1)
                    connec = obj.localGlobalDofConnec{jdom,idom};
                    Res = obj.local2global(RHS{jdom,idom},connec);
                    Rg = Res-rCoarse;
                    Rbal{jdom,idom} = obj.global2local(Rg,connec);
                end
            end          
        end

%          function Rbal = balanceResidual(obj,u,RHS)
%             R=zeros(size(obj.interfaceDof,1),1);
% %             R2=zeros(size(obj.interfaceDof,1),1);
%             unodal = reshape(u',1,[])';
%             nint = size(obj.interfaceDof,3);
%             for iint = 1:nint
%                 ndom = size(obj.interfaceDof(:,:,iint),2);
%                 for idom = 1:ndom
%                     dom = obj.interfaceDom(iint,idom);
%                     row = ceil(dom/obj.nSubdomains(1));
%                     col = dom-(row-1)*obj.nSubdomains(1);
% %                     unodal = reshape(u{row,col}',1,[])';
%                     dof = obj.interfaceDof(:,idom,iint);
% %                      R = R + w(idom)*(RHS{row,col}(dof)-obj.LHS{row,col}(dof,:)*unodal ...
% %                     + obj.LHS{row,col}(dof,dof)*unodal(dof));
%                     connec = obj.localGlobalDofConnec{row,col};
%                     Res = obj.local2global(RHS{row,col},connec);
%                     Rg = Res-obj.domainLHS*unodal;
%                     Rl = obj.global2local(Rg,connec);
%                     Rbal(:,idom,iint)=Rl(dof);
%                 end
% %                 Rbal(:,,iint)=R;
%             end
%         end

        function uInt = computeInterfaceDisp(obj,u)
            nint = size(obj.interfaceDof,3);
            uInt = zeros(size(obj.interfaceDof,1),nint);
            w = [obj.weight,1- obj.weight];
            for iint = 1:nint
                ndom = size(obj.interfaceDof(:,:,iint),2);
                for idom = 1:ndom
                    dom = obj.interfaceDom(iint,idom);
                    row = ceil(dom/obj.nSubdomains(1));
                    col = dom-(row-1)*obj.nSubdomains(1);
                    unodal = reshape(u{row,col}',1,[])';
                    dof = obj.interfaceDof(:,idom,iint);
                    uInt(:,iint) = uInt(:,iint) + w(idom)*unodal(dof);
                end
            end
%             w = [obj.weight,1- obj.weight];
%             %             w  = [1,1];
%             uInt=zeros(size(obj.interfaceDof,1),1);
%             for idom=1:obj.nSubdomains(1)
%                 unodal = reshape(u{idom}',1,[])';
%                 interfaceDOF = obj.interfaceDof(:,idom);
%                 uInt = uInt + w(idom)*unodal(interfaceDOF);
%             end
        end

        function uIntNew = updateDirichletValues(obj,u)
            nint = size(obj.interfaceDof,3);
            for iint = 1:nint
                ndom = size(obj.interfaceDof(:,:,iint),2);
                for idom = 1:ndom
                    dom = obj.interfaceDom(iint,idom);
                    row = ceil(dom/obj.nSubdomains(1));
                    col = dom-(row-1)*obj.nSubdomains(1);
%                     ind = obj.interfaceDirichletDof(:,idom,iint);
                    dof = obj.interfaceDof(:,idom,iint);
                    bD{row,col} = obj.boundaryConditions.dirichletStep{row,col};
%                     bD{row,col}.dirichlet_values(ind) = bD{row,col}.dirichlet_values(ind)+ obj.theta*u(:,iint);
%                     bc  = obj.boundaryConditions.dirichletStep{row,col};
                    values = reshape(bD{row,col}.dirichletFun.fValues',1,[])';
                    values(dof) = values(dof) + obj.theta*u(:,iint);
                    values = reshape(values,obj.ndimf,[])';
                    bD{row,col}.dirichletFun.fValues = values;
                end
                ui = reshape(bD{row,col}.dirichletFun.fValues',1,[])';
                uIntNew(:,iint) = ui(dof);
            end
            obj.boundaryConditions.dirichletStep = bD;
%             uIntNew = bD{end}.dirichlet_values(obj.interfaceDirichletDof(:,end));
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
            s.ndimf = 3;
            xF = P1Function(s);
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
            nbasis = size(basis,2);
            for ibasis=1:nbasis
                    nfields = size(basis{ibasis},2);
                    for ifields = 1:nfields
                        values = basis{ibasis}(:,ifields);
                        s.fValues = reshape(values,2,[])';
                        s.mesh = obj.meshDomain;
                        s.fValues(:,end+1) = 0;
                        s.ndimf = 3;
                        fun{ifields} = P1Function(s);
                    end
                    %                 a.fun=fun;
                    mesh=obj.meshDomain;
                    filename = ['domain',num2str(ibasis)];
                    obj.print(filename,software,funNames,fun,mesh)
            end
        end

        function bc = bcDirichletStepLibrary(obj)
            isLeft   = @(coor) (abs(coor(:,1) - min(coor(:,1)))   < 1e-12);
            isRight  = @(coor) (abs(coor(:,1) - max(coor(:,1)))   < 1e-12);
            isBottom = @(coor) (abs(coor(:,2) - min(coor(:,2)))   < 1e-12);
            isTop    = @(coor) (abs(coor(:,2) - max(coor(:,2)))   < 1e-12);

            bc.dirichletBottomLeft.domain    = @(coor) isLeft(coor)|isRight(coor) ;
            bc.dirichletBottomLeft.direction = [1,2];
            bc.dirichletBottomLeft.value     = 0;

            bc.dirichletBottomMiddle.domain    = @(coor) isLeft(coor) | isRight(coor) ;
            bc.dirichletBottomMiddle.direction = [1,2];
            bc.dirichletBottomMiddle.value     = 0;

            bc.dirichletBottomRight.domain    = @(coor) isLeft(coor)  ;
            bc.dirichletBottomRight.direction = [1,2];
            bc.dirichletBottomRight.value     = 0;

% 
%             bc.dirichletBottomLeft.boundaryId=[1,2];
%             bc.dirichletBottomLeft.dof{1}=[1,2];
%             bc.dirichletBottomLeft.dof{2}=[1,2];
%             bc.dirichletBottomLeft.value{1}=[0,0];
%             bc.dirichletBottomLeft.value{2}=[0,0];
% %             
% %             bc.dirichletBottomLeft.boundaryId=[1,2,4];
% %             bc.dirichletBottomLeft.dof{1}=[1,2];
% %             bc.dirichletBottomLeft.dof{2}=[1,2];
% %             bc.dirichletBottomLeft.dof{3}=[1,2];
% %             bc.dirichletBottomLeft.value{1}=[0,0];
% %             bc.dirichletBottomLeft.value{2}=[0,0];
% %             bc.dirichletBottomLeft.value{3}=[0,0];
%             %             bc.newmanNo=[];
% 
%             bc.dirichletBottomMiddle.boundaryId=[1,2];
%             bc.dirichletBottomMiddle.dof{1}=[1,2];
%             bc.dirichletBottomMiddle.dof{2}=[1,2];
%             bc.dirichletBottomMiddle.value{1}=[0,0];
%             bc.dirichletBottomMiddle.value{2}=[0,0];
% 
% %             bc.dirichletBottomMiddle.boundaryId=[1,2,4];
% %             bc.dirichletBottomMiddle.dof{1}=[1,2];
% %             bc.dirichletBottomMiddle.dof{2}=[1,2];
% %             bc.dirichletBottomMiddle.dof{3}=[1,2];
% %             bc.dirichletBottomMiddle.value{1}=[0,0];
% %             bc.dirichletBottomMiddle.value{2}=[0,0];
% %             bc.dirichletBottomMiddle.value{3}=[0,0];
% 
%             bc.dirichletBottomRight.boundaryId=[1];
%             bc.dirichletBottomRight.dof{1}=[1,2];
% %             bc.dirichletBottomRight.dof{2}=[1,2];
%             bc.dirichletBottomRight.value{1}=[0,0];
% %             bc.dirichletBottomRight.value{2}=[0,0];
% 
% %             bc.dirichletBottomRight.boundaryId=[1,4];
% %             bc.dirichletBottomRight.dof{1}=[1,2];
% %             bc.dirichletBottomRight.dof{2}=[1,2];
% %             bc.dirichletBottomRight.value{1}=[0,0];
% %             bc.dirichletBottomRight.value{2}=[0,0];

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

            bc.dirichletLeftMiddle.boundaryId=[1,2,3,4];
            bc.dirichletLeftMiddle.dof{1}=[1,2];
            bc.dirichletLeftMiddle.dof{2}=[1,2];
            bc.dirichletLeftMiddle.dof{3}=[1,2];
            bc.dirichletLeftMiddle.dof{4}=[1,2];
            bc.dirichletLeftMiddle.value{1}=[0,0];
            bc.dirichletLeftMiddle.value{2}=[0,0];
            bc.dirichletLeftMiddle.value{3}=[0,0];
            bc.dirichletLeftMiddle.value{4}=[0,0];

            bc.dirichletRightMiddle.boundaryId=[1,3,4];
            bc.dirichletRightMiddle.dof{1}=[1,2];
            bc.dirichletRightMiddle.dof{2}=[1,2];
            bc.dirichletRightMiddle.dof{3}=[1,2];
            bc.dirichletRightMiddle.value{1}=[0,0];
            bc.dirichletRightMiddle.value{2}=[0,0];
            bc.dirichletRightMiddle.value{3}=[0,0];

            bc.dirichletMiddleMiddle.boundaryId=[1,2,3,4];
            bc.dirichletMiddleMiddle.dof{1}=[1,2];
            bc.dirichletMiddleMiddle.dof{2}=[1,2];
            bc.dirichletMiddleMiddle.dof{3}=[1,2];
            bc.dirichletMiddleMiddle.dof{4}=[1,2];
            bc.dirichletMiddleMiddle.value{1}=[0,0];
            bc.dirichletMiddleMiddle.value{2}=[0,0];
            bc.dirichletMiddleMiddle.value{3}=[0,0];
            bc.dirichletMiddleMiddle.value{4}=[0,0];
        end

        function bc = bcNeumannStepLibrary(obj)
            isLeft   = @(coor) (abs(coor(:,1) - min(coor(:,1)))   < 1e-12);
            isRight  = @(coor) (abs(coor(:,1) - max(coor(:,1)))   < 1e-12);
            isBottom = @(coor) (abs(coor(:,2) - min(coor(:,2)))   < 1e-12);
            isTop    = @(coor) (abs(coor(:,2) - max(coor(:,2)))   < 1e-12);

            bc.neumannBottomLeft.domain    = @(coor) isRight(coor) ;
            bc.neumannBottomLeft.direction = [1,2];
            bc.neumannBottomLeft.value     = 0;

            bc.neumannBottomMiddle.domain    = @(coor) isLeft(coor) | isRight(coor) ;
            bc.neumannBottomMiddle.direction = [1,2];
            bc.neumannBottomMiddle.value     = 0;

            bc.neumannBottomRight.domain    = @(coor) isLeft(coor) ;
            bc.neumannBottomRight.direction = [1,2];
            bc.neumannBottomRight.value     = 0;

% 
%             bc.neumannBottomLeft.boundaryId=[2];
%             bc.neumannBottomLeft.dof{1}=[1,2];
%             bc.neumannBottomLeft.value{1}=[0,0];
% % %             %             bc.newmanNo=[];
%             
% %             bc.neumannBottomLeft.boundaryId=[2,4];
% %             bc.neumannBottomLeft.dof{1}=[1,2];
% %             bc.neumannBottomLeft.dof{2}=[1,2];
% %             bc.neumannBottomLeft.value{1}=[0,0];
% %             bc.neumannBottomLeft.value{2}=[0,0];
%             %             bc.newmanNo=[];
% 
%             bc.neumannBottomMiddle.boundaryId=[1,2];
%             bc.neumannBottomMiddle.dof{1}=[1,2];
%             bc.neumannBottomMiddle.dof{2}=[1,2];
%             bc.neumannBottomMiddle.value{1}=[0,0];
%             bc.neumannBottomMiddle.value{2}=[0,0];
%             
% %             bc.neumannBottomMiddle.boundaryId=[1,2,4];
% %             bc.neumannBottomMiddle.dof{1}=[1,2];
% %             bc.neumannBottomMiddle.dof{2}=[1,2];
% %             bc.neumannBottomMiddle.dof{3}=[1,2];
% %             bc.neumannBottomMiddle.value{1}=[0,0];
% %             bc.neumannBottomMiddle.value{2}=[0,0];
% %             bc.neumannBottomMiddle.value{3}=[0,0];
% 
%             bc.neumannBottomRight.boundaryId=[1];
%             bc.neumannBottomRight.dof{1}=[1,2];
%             bc.neumannBottomRight.value{1}=[0,0];
% 
% %             bc.neumannBottomRight.boundaryId=[1,4];
% %             bc.neumannBottomRight.dof{1}=[1,2];
% %             bc.neumannBottomRight.dof{2}=[1,2];
% %             bc.neumannBottomRight.value{1}=[0,0];
% %             bc.neumannBottomRight.value{2}=[0,0];

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

            bc.neumannLeftMiddle.boundaryId=[2,3,4];
            bc.neumannLeftMiddle.dof{1}=[1,2];
            bc.neumannLeftMiddle.dof{2}=[1,2];
            bc.neumannLeftMiddle.dof{3}=[1,2];
            bc.neumannLeftMiddle.value{1}=[0,0];
            bc.neumannLeftMiddle.value{2}=[0,0];
            bc.neumannLeftMiddle.value{3}=[0,0];

            bc.neumannRightMiddle.boundaryId=[1,3,4];
            bc.neumannRightMiddle.dof{1}=[1,2];
            bc.neumannRightMiddle.dof{2}=[1,2];
            bc.neumannRightMiddle.dof{3}=[1,2];
            bc.neumannRightMiddle.value{1}=[0,0];
            bc.neumannRightMiddle.value{2}=[0,0];
            bc.neumannRightMiddle.value{3}=[0,0];

            bc.neumannMiddleMiddle.boundaryId=[1,2,3,4];
            bc.neumannMiddleMiddle.dof{1}=[1,2];
            bc.neumannMiddleMiddle.dof{2}=[1,2];
            bc.neumannMiddleMiddle.dof{3}=[1,2];
            bc.neumannMiddleMiddle.dof{4}=[1,2];
            bc.neumannMiddleMiddle.value{1}=[0,0];
            bc.neumannMiddleMiddle.value{2}=[0,0];
            bc.neumannMiddleMiddle.value{3}=[0,0];
            bc.neumannMiddleMiddle.value{4}=[0,0];
        end
    end

    
end
