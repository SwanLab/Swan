classdef EIFEMtesting < handle

    properties (Access = public)

    end

    properties (Access = private)

    end

    properties (Access = private)
        EIFEMfilename
        meshDomain
        boundaryConditions
        bcApplier
        solverCase
        material
        meshReference
        interfaceMeshReference
        meshSubDomain
        nSubdomains
        ninterfaces
        interfaceMeshSubDomain
        globalMeshConnecSubDomain
        interfaceMeshConnecSubDomain
        subDomainContact
        cornerNodes
        quad
        interfaceConnec
        locGlobConnec
        localGlobalDofConnec

        displacementFun
        LHS
        RHS
        scale
        ndimf    
        functionType
        EIFEM
    end

    methods (Access = public)

        function obj = EIFEMtesting()
            close all
            obj.init();
            obj.createReferenceMesh();
           
            s.nsubdomains   = obj.nSubdomains; %nx ny
            s.meshReference = obj.meshReference;         
            m = MeshCreatorFromRVE(s);
            [obj.meshDomain,obj.meshSubDomain,obj.interfaceConnec,~,obj.locGlobConnec] = m.create();

            obj.displacementFun      = LagrangianFunction.create(obj.meshDomain, obj.ndimf,obj.functionType);
            [obj.boundaryConditions,Dir,PL] = obj.createBoundaryConditions(obj.meshDomain);
            ss.mesh                  = obj.meshDomain;
            ss.boundaryConditions    = obj.boundaryConditions;
            obj.bcApplier            = BCApplier(ss);
            obj.localGlobalDofConnec = obj.createlocalGlobalDofConnec();
%             obj.quad               = Quadrature.set(obj.meshDomain.type);
%             obj.quad.computeQuadrature('QUADRATIC');
%             obj.createDomainMaterial();
%               obj.computeForces();
    
            obj.material = obj.createMaterial(obj.meshDomain);
            obj.LHS      = obj.computeStiffnessMatrix();
            obj.RHS      = obj.computeForces(obj.LHS);
           
            cMesh           = createCoarseMesh(obj);
            s.nsubdomains   = obj.nSubdomains; %nx ny
            s.meshReference = cMesh;         
            mRVECoarse      = MeshCreatorFromRVE(s);
            [meshDomainCoarse,meshSubDomainCoarse,interfaceConnecCoarse] = mRVECoarse.create();

            obj.EIFEM = obj.createEIFEM(meshDomainCoarse,Dir);

            LHS = obj.bcApplier.fullToReducedMatrixDirichlet(obj.LHS);
            RHS = obj.bcApplier.fullToReducedVectorDirichlet(obj.RHS);

            u = obj.solver(LHS,RHS);

% %             obj.obtainCornerNodes();
% %             fineMesh = MeshFromRVE
%             obj.createSubDomainMeshes();
%             obj.createInterfaceSubDomainMeshes();
%             obj.createDomainMesh();

            
%             s.referenceMesh = obj.referenceMesh;
%             mC = MeshCreatorFromSubmeshes();
%             obj.meshDomain = mC.mesh;

%             preconditioner = obj.createPreconditioner(mC.submeshes);

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
            obj.nSubdomains  = [2 1]; %nx ny
            obj.scale        = 'MACRO';
            obj.ndimf        = 2;
            obj.functionType = 'P1';
            obj.solverCase   = 'REDUCED';
            obj.EIFEMfilename  = '/home/raul/Documents/Thesis/EIFEM/EXAMPLE/EIFE_LIBRARY/DEF_Q4por_1.mat';
        end

        function createReferenceMesh(obj)
%             filename   = 'lattice_ex1';
%             a.fileName = obj.EIFEMfilename;
%             femD       = FemDataContainer(a);
%             mS         = femD.mesh;
%             bS         = mS.createBoundaryMesh();
%              % Generate coordinates
%             x1 = linspace(0,1,5);
%             x2 = linspace(0,1,5);
%             % Create the grid
%             [xv,yv] = meshgrid(x1,x2);
%             % Triangulate the mesh to obtain coordinates and connectivities
%             [F,coord] = mesh2tri(xv,yv,zeros(size(xv)),'x');
% 
%             s.coord    = coord(:,1:2);
%             s.connec   = F;
%             mS         = Mesh.create(s);
%             bS         = mS.createBoundaryMesh();
            load(obj.EIFEMfilename);
            s.coord    = EIFEoper.MESH.COOR;
            s.connec   = EIFEoper.MESH.CN;
            mS         = Mesh.create(s);
            bS         = mS.createBoundaryMesh();
            
            obj.meshReference = mS;
            obj.interfaceMeshReference = bS;
            obj.ninterfaces = length(bS);
        end

        function cMesh = createCoarseMesh(obj)
            xmax = max(obj.meshReference.coord(:,1));
            xmin = min(obj.meshReference.coord(:,1));
            ymax = max(obj.meshReference.coord(:,2));
            ymin = min(obj.meshReference.coord(:,2));
            coord(1,1) = xmax;
            coord(1,2) = ymin;
            coord(2,1) = xmax;
            coord(2,2) = ymax;
            coord(3,1) = xmin;
            coord(3,2) = ymax;
            coord(4,1) = xmin;
            coord(4,2) = ymin;
            connec = [1 2 3 4];
            s.coord = coord;
            s.connec = connec;
            cMesh = Mesh.create(s);
        end

        function createDomainMaterial(obj)
%             ngaus = 1;        
            m = obj.meshDomain;                      
            obj.material = obj.createMaterial(m);
        end
        
        function L = computeReferenceMeshLength(obj)
            coord = obj.meshReference.coord;
            Lx = max(coord(:,1));
            Ly = max(coord(:,2));
            L = [Lx Ly];
        end     

        function [Dir,PL] = createRawBoundaryConditions(obj)
            isLeft   = @(coor) (abs(coor(:,1) - min(coor(:,1)))   < 1e-12);
            isRight  = @(coor) (abs(coor(:,1) - max(coor(:,1)))   < 1e-12);
            %             isMiddle = @(coor) (abs(coor(:,2) - max(coor(:,2)/2)) == 0);
            Dir.domain    = @(coor) isLeft(coor);
            Dir.direction = [1,2];
            Dir.value     = 0;

            PL.domain    = @(coor) isRight(coor);
            PL.direction = 2;
            PL.value     = -1;
        end        

         function [bc,Dir,PL] = createBoundaryConditions(obj,mesh)
            [Dir,PL]  = obj.createRawBoundaryConditions();
            dirichlet = DirichletCondition(mesh,Dir);
            pointload = PointLoad(mesh,PL);
             % need this because force applied in the face not in a point
            pointload.values        = pointload.values/size(pointload.dofs,1);
            fvalues                 = zeros(mesh.nnodes*obj.ndimf,1);
            fvalues(pointload.dofs) = pointload.values;
            fvalues                 = reshape(fvalues,obj.ndimf,[])';
            pointload.fun.fValues   = fvalues;

            s.pointloadFun = pointload;
            s.dirichletFun = dirichlet;
            s.periodicFun  =[];
            s.mesh         = mesh;
            bc             = BoundaryConditions(s);
        end

        function [dirichlet,pointload] = createBc(obj,boundaryMesh,dirchletBc,newmanBc)
            dirichlet = obj.createBondaryCondition(boundaryMesh,dirchletBc);
            pointload = obj.createBondaryCondition(boundaryMesh,newmanBc);
        end

        function cond = createBondaryCondition(obj,bM,condition)
            nbound = length(condition.boundaryId);
            cond = zeros(1,3);
            for ibound=1:nbound
                ncond  = length(condition.dof(nbound,:));
                nodeId= unique(bM{condition.boundaryId(ibound)}.globalConnec);
                nbd   = length(nodeId);
                condition.value{ibound} = condition.value{ibound}/nbd;
                for icond=1:ncond
                    bdcond= [nodeId, repmat(condition.dof(icond),[nbd,1]), repmat(condition.value(icond),[nbd,1])];
                    cond=[cond;bdcond];
                end
            end
            cond = cond(2:end,:);
        end

%         function material = createMaterial(obj,mesh)
%             I = ones(mesh.nelem,obj.quad.ngaus);
%             s.ptype = 'ELASTIC';
%             s.pdim  = '2D';
%             s.nelem = mesh.nelem;
%             s.mesh  = mesh;
%             s.kappa = .9107*I;
%             s.mu    = .3446*I;
%             mat = Material.create(s);
%             mat.compute(s);
%             material = mat;
%         end
        
        function [young,poisson] = computeElasticProperties(obj,mesh)
            E1  = 1;
            nu1 = 1/3;
            E   = AnalyticalFunction.create(@(x) E1*ones(size(squeeze(x(1,:,:)))),1,mesh);
            nu  = AnalyticalFunction.create(@(x) nu1*ones(size(squeeze(x(1,:,:)))),1,mesh);
            young   = E;
            poisson = nu;
        end

        function material = createMaterial(obj,mesh)
            [young,poisson] = obj.computeElasticProperties(mesh);
            s.type    = 'ISOTROPIC';
            s.ptype   = 'ELASTIC';
            s.ndim    = mesh.ndim;
            s.young   = young;
            s.poisson = poisson;
            tensor    = Material.create(s);
            material  = tensor;
        end    

        function dim = getFunDims(obj)
            d.ndimf  = obj.displacementFun.ndimf;
            d.nnodes = size(obj.displacementFun.fValues, 1);
            d.ndofs  = d.nnodes*d.ndimf;
            d.nnodeElem = obj.meshDomain.nnodeElem; % should come from interp..
            d.ndofsElem = d.nnodeElem*d.ndimf;
            dim = d;
        end

%         function computeStiffnessMatrix(obj)
%             s.type     = 'ElasticStiffnessMatrix';
%             s.mesh     = obj.meshDomain;
%             s.fun      = obj.displacementFun;
%             s.material = obj.material;
%             lhs = LHSintegrator.create(s);
%             obj.LHS = lhs.compute();
%         end

        function LHS = computeStiffnessMatrix(obj)
            s.type     = 'ElasticStiffnessMatrix';
            s.mesh     = obj.meshDomain;
            s.test     = LagrangianFunction.create(s.mesh,obj.ndimf, obj.functionType);
            s.trial    = obj.displacementFun;
            s.material = obj.material;
            s.quadratureOrder = 2;
            lhs = LHSintegrator.create(s);
            LHS = lhs.compute();
        end

%         function computeForces(obj)
%             s.type = 'Elastic';
%             s.scale    = obj.scale;
%             s.dim      = obj.getFunDims();
%             s.BC       = obj.boundaryConditions;
%             s.mesh     = obj.meshDomain;
%             s.material = obj.material;
% %             s.globalConnec = obj.displacementField.connec;
%             s.globalConnec = obj.meshDomain.connec;
%             RHSint = RHSintegrator.create(s);
%             rhs = RHSint.compute();
%             R = RHSint.computeReactions(obj.LHS);
% %             obj.variables.fext = rhs + R;
%             obj.RHS = rhs;
%         end

        function forces = computeForces(obj,stiffness)
            s.type      = 'Elastic';
            s.scale     = 'MACRO';
%             s.dim       = obj.getFunDims();
            s.dim.ndofs = obj.displacementFun.nDofs;
            s.BC        = obj.boundaryConditions;
            s.mesh      = obj.meshDomain;
            s.material  = obj.material;
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
         
        function Gvec = local2global(obj,Lvec,locGlobConnec)
            ndimf  = obj.displacementFun.ndimf;
            Gvec   = zeros(obj.meshDomain.nnodes*ndimf,1);
            Gvec(locGlobConnec(:,1)) = Lvec(locGlobConnec(:,2));
        end

        function Lvec = global2local(obj,Gvec)
            ndimf  = obj.displacementFun.ndimf;
            Lvec   = zeros(obj.meshReference.nnodes*ndimf,obj.nSubdomains(1)*obj.nSubdomains(2));
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
            RVE         = TrainedRVE(obj.EIFEMfilename);
            s1.RVE      = RVE;
            s1.mesh     = meshDomainCoarse;
            s1.DirCond  = Dir;
            eifem       = EIFEM(s1);
        end

        function u = solver(obj,LHS,RHS)
            tol=1e-8;
            e=1;
            iter=1;
            u = zeros(length(RHS),1);
            while e(iter)>tol
                R = RHS - LHS*u;
                RG = obj.bcApplier.reducedToFullVectorDirichlet(R);
                RGsbd = obj.global2local(RG);


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

    end
end
