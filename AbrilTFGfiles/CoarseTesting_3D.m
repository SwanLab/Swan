classdef CoarseTesting_3D< handle

    properties (Access = public)
        SolExact
        Sol
        residualPCG
        errPCG
        errAnormPCG
        residualCG
        errCG
        errAnormCG
        residualILU
    end

    properties (Access = private)
        nSubdomains
        r
        centroids
        ic
        icr
        lg
        bs
        meshDomain
        referenceMesh
        subdomainMeshes
        cellMesh
        discMesh
        boundaryConditions
        bcApplier
        ddDofManager

        tolSameNode
        data
        params

        xmin 
        xmax 
        ymin 
        ymax  
        zmin
        zmax
        designVariable
        materialInterpolator
        unfittedMesh
        fileNameEIFEM

    end


    methods (Access = public)

        function obj = CoarseTesting_3D(cParams)
            obj.init(cParams)
        end

        function compute(obj)
            mR  = obj.createReferenceMesh();
            obj.referenceMesh=mR;
            bS  = mR.createBoundaryMesh();
            obj.createMeshDomain(mR);
            [bC,dir] = obj.createBoundaryConditions(obj.meshDomain);
            obj.boundaryConditions = bC;
            obj.createBCapplier()
            [LHS,RHS,~] = obj.createElasticProblem();

            %EXACT SOLUTION
            tic
            LHSf   = @(x) LHS*x;
            RHSf   = RHS;
            Usol   = LHS\RHS;
            Ufull  = obj.bcApplier.reducedToFullVectorDirichlet(Usol); 
            t_direct=toc

            % PRECONDITIONERS
            % Milu        = obj.createILUpreconditioner(LHS);
            Mid         = @(r)r;
            switch obj.params.Option
                case {'Dataset','NN','Hybrid'}
                    Mcoarse     = obj.createCoarseNNPreconditioner(mR,dir,obj.ic,obj.lg,bS,obj.icr,obj.discMesh);
                    Mmult        = @(r) Preconditioner.multiplePrec(r,LHSf,Milu,Mcoarse,Milu);
                case {'DIRECT','Direct'}
                    Meifem       = obj.createEIFEMPreconditioner(dir,obj.ic,obj.lg,bS,obj.icr,obj.discMesh);
                    Mmult        = @(r) Preconditioner.multiplePrec(r,LHSf,Milu,Meifem,Milu);
            end

            tol = 1e-8;
            x0  = zeros(size(RHSf));

            % tic %SOLVE THE CASE WITH STANDARD CG
            % [~,obj.residualCG,errCG, errCG] = PCG.solve(LHSf,RHSf,x0,Mid,tol,Usol,obj.meshDomain,obj.bcApplier);
            % t_CG=toc
            % tic  % SOLVE THE CASE WITH CG+ ILU
            % [~,obj.residualILU,errILU, errAnormILU] = PCG.solve(LHSf,RHSf,x0,Milu,tol,Usol,obj.meshDomain,obj.bcApplier);
            % t_ILU=toc
           % SOLVE THE CASE WITH PRECONDITIONING ILU+EIFEM+ILU
            tic
            [uPCG,obj.residualPCG,obj.errPCG,obj.errAnormPCG] = PCG.solve(LHSf,RHSf,x0,Mmult,tol,Usol,obj.meshDomain,obj.bcApplier);
            t_PCG=toc
            xFull = obj.bcApplier.reducedToFullVectorDirichlet(uPCG);
            

            uDomain = obj.bcApplier.reducedToFullVectorDirichlet(uPCG);
            % uDomain = obj.ddDofManager.global2local(uDomain);
            
            % % LAGRANGIAN FUN SOLUTIONS
            % s.mesh     = obj.meshDomain;
            % s.ndimf    = obj.meshDomain.ndim;
            % s.order    = 'P1';
            % s.fValues  = reshape(xFull,3,[])';
            % obj.Sol    = LagrangianFunction(s); %Preconditioned sol  
            % s.fValues = reshape(Ufull,3[])';
            % obj.SolExact=LagrangianFunction(s); %Exact sol

            % obj.print(uPCG,"SolPCG");
            %CoarsePlotSolution(uFun, obj.meshDomain, obj.bcApplier,'TestCoarseAbril', obj.r, obj.centroids);
            %CoarsePlotSolution(RealFun, obj.meshDomain, obj.bcApplier,'TestRealAbril', obj.r, obj.centroids);

        end

        function PlotSolution(obj)
            figure('Position',[65 280 1800 600])
            tiledlayout(1,3)

            nexttile
            plot(obj.residualPCG,'linewidth',2)
            hold on
            plot(obj.residualILU,'linewidth',2)
            set(gca, 'YScale', 'log')
            xlabel('Iteration')
            ylabel('Residual')
            title("Residual")
            legend({'PCG','CG'})

            nexttile
            plot(obj.errPCG,'linewidth',2)
            hold on
            plot(errILU,'linewidth',2)
            set(gca, 'YScale', 'log')
            xlabel('Iteration')
            ylabel('||error||_{L2}')
            title("error")
            legend({'PCG','CG'})

            nexttile
            plot(obj.errAnormPCG,'linewidth',2)
            hold on
            plot(errAnormILU,'linewidth',2)
            set(gca, 'YScale', 'log')
            xlabel('Iteration')
            ylabel('Energy norm')
            title("Err Anorm")
            legend({'PCG','CG'})
        end

        function print(obj,sol,fileName)
                z.mesh      = obj.meshDomain;
                z.order     = 'P1';
                z.fValues   = reshape(sol,z.mesh.ndim,[])';
                uFeFun = LagrangianFunction(z);%
                uMeshFun = obj.unfittedMesh.obtainFunctionAtUnfittedMesh(uFeFun);
                
                fvalues = [uMeshFun.innerMeshFunction.fValues;
                    uMeshFun.innerCutMeshFunction.fValues];

                s.coord = [uMeshFun.innerMeshFunction.mesh.coord;
                    uMeshFun.innerCutMeshFunction.mesh.coord];

                s.connec = [uMeshFun.innerMeshFunction.mesh.connec;
                    uMeshFun.innerCutMeshFunction.mesh.connec  + max(uMeshFun.innerMeshFunction.mesh.connec(:))];
                
                mh = Mesh.create(s);
                ss.mesh = mh;
                ss.fValues = fvalues;
                ss.order = 'P1';
                ss.ndimf = size(fvalues,2)
                u = LagrangianFunction(ss);
                u.print(fileName,'Paraview')
        end

    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods (Access = private)

        function init(obj,cParams)
            p.Training      = cParams.Training;    % 'EIFEM'/'Multiscale'
            p.Sampling      = cParams.Sampling;    % 'Isolated'/'Oversampling'
            p.Inclusion     = cParams.Inclusion;   % 'Hole'/'Material'/'HoleRaul'   --> Hole: just for constant r
            p.Option        = cParams.Option;      % 'Dataset'/'NN'/'Direct'
            p.nelem         = cParams.nelem;       %  Mesh refining
            obj.params      = p;
            obj.r           = cParams.r;
            [Ny,Nx,Nz] = size(obj.r);
            obj.nSubdomains = [Nx Ny Nz];
            % obj.nSubdomains = [35 1 1];    % UNCOMMENT JUST FOR AIRFOIL
            obj.tolSameNode = 1e-10;   % 1E-10--> general case   1E-6 --> airfoil
            obj.fileNameEIFEM = cParams.fileNameEIFEM;
        end


        function createMeshDomain(obj,mR)
            s.nsubdomains   = obj.nSubdomains; %  num along x and y
            s.meshReference = mR;
            s.tolSameNode = obj.tolSameNode;
            m = MeshCreatorFromRVE3D(s);
            [mD,mSb,iC,~,lG,iCR,discmesh] = m.create();
            obj.meshDomain      = mD;        
            obj.subdomainMeshes = mSb;  
            obj.ic              = iC;   
            obj.icr             = iCR;  
            obj.lg              = lG;   
            obj.discMesh        = discmesh;
        end

        function mS = createReferenceMesh(obj)
            mS = obj.createStructuredMesh();
            % mS = obj.importGIDMesh();
        end


        function mS = createStructuredMesh(obj)
            n =obj.params.nelem;
            m = TetraMesh(1,1,1,n,n,n);

            s.coord=m.coord;
            s.connec=m.connec;

            maxC= max(s.coord);
            minC = min(s.coord);

            obj.xmin = minC(1);
            obj.xmax = maxC(1);
            obj.ymin = minC(2);
            obj.ymax = maxC(2);
            obj.zmin = minC(3);
            obj.zmax = maxC(3);

            % delta = 1e-6;
            % tol   = 1e-12;
            % % Detect boundaries per direction
            % onMin = abs(s.coord - minC) < tol;
            % onMax = abs(s.coord - maxC) < tol;
            % 
            % isBoundary = onMin | onMax;
            % 
            % % Count how many boundary planes each node lies on
            % numBoundaries = sum(isBoundary, 2);
            % 
            % % Edge nodes = exactly 2 boundaries
            % isEdge = numBoundaries == 2;
            % 
            % % Displacement field
            % dispVec = zeros(size(s.coord));
            % dispVec(onMax) = -delta;
            % dispVec(onMin) =  delta;
            % 
            % % Apply only to edges
            % s.coord(isEdge,:) = s.coord(isEdge,:) + dispVec(isEdge,:);
            % 
            % isEdgeOrCorner = numBoundaries > 2;
            % s.coord(isEdgeOrCorner,:) = s.coord(isEdgeOrCorner,:) + dispVec(isEdgeOrCorner,:);

            delta=1E-8;

            s.coord(s.coord(:,1)== maxC(1) & s.coord(:,3)==maxC(3),:) =...
                s.coord(s.coord(:,1)== maxC(1) & s.coord(:,3)==maxC(3),:)-[0,0,delta];

            s.coord(s.coord(:,1)== maxC(1) & s.coord(:,3)==minC(3),:) =...
                s.coord(s.coord(:,1)== maxC(1) & s.coord(:,3)==minC(3),:)+[0,0,delta];

            s.coord(s.coord(:,1)== minC(1) & s.coord(:,3)==maxC(3),:) =...
                s.coord(s.coord(:,1)== minC(1) & s.coord(:,3)==maxC(3),:)-[0,0,delta];

            s.coord(s.coord(:,1)== minC(1) & s.coord(:,3)==minC(3),:) =...
                s.coord(s.coord(:,1)== minC(1) & s.coord(:,3)==minC(3),:)+[0,0,delta];

            s.coord(s.coord(:,2)== maxC(2) & s.coord(:,3)==maxC(3),:) =...
                s.coord(s.coord(:,2)== maxC(2) & s.coord(:,3)==maxC(3),:)-[0,0,delta];

            s.coord(s.coord(:,2)== maxC(2) & s.coord(:,3)==minC(3),:) =...
                s.coord(s.coord(:,2)== maxC(2) & s.coord(:,3)==minC(3),:)+[0,0,delta];

            s.coord(s.coord(:,2)== minC(2) & s.coord(:,3)==maxC(3),:) =...
                s.coord(s.coord(:,2)== minC(2) & s.coord(:,3)==maxC(3),:)-[0,0,delta];

            s.coord(s.coord(:,2)== minC(2) & s.coord(:,3)==minC(3),:) =...
                s.coord(s.coord(:,2)== minC(2) & s.coord(:,3)==minC(3),:)+[0,0,delta];

            s.coord(s.coord(:,1)== maxC(1) & s.coord(:,2)==maxC(2),:) =...
                s.coord(s.coord(:,1)== maxC(1) & s.coord(:,2)==maxC(2),:)-[0,delta,0];

            s.coord(s.coord(:,1)== maxC(1) & s.coord(:,2)==minC(2),:) =...
                s.coord(s.coord(:,1)== maxC(1) & s.coord(:,2)==minC(2),:)+[0,delta,0];

            s.coord(s.coord(:,1)== minC(1) & s.coord(:,2)==maxC(2),:) =...
                s.coord(s.coord(:,1)== minC(1) & s.coord(:,2)==maxC(2),:)-[0,delta,0];

            s.coord(s.coord(:,1)== minC(1) & s.coord(:,2)==minC(2),:) =...
                s.coord(s.coord(:,1)== minC(1) & s.coord(:,2)==minC(2),:)+[0,delta,0];

            mS = Mesh.create(s);

        end


        function mS=importGIDMesh(obj)
            filename = 'DEF_Q8_wing_1.mat';
            load(filename);
            s.coord    = EIFEoper.MESH.COOR;
            s.connec   = EIFEoper.MESH.CN;

            maxC= max(s.coord);
            minC = min(s.coord);

            obj.xmin = minC(1);
            obj.xmax = maxC(1);
            obj.ymin = minC(2);
            obj.ymax = maxC(2);
            obj.zmin = minC(3);
            obj.zmax = maxC(3);

            delta=1E-5;  % 1E-9--> general case   1E-5 --> airfoil

            s.coord(s.coord(:,1)== maxC(1) & s.coord(:,3)==maxC(3),:) =...
                s.coord(s.coord(:,1)== maxC(1) & s.coord(:,3)==maxC(3),:)-[0,0,delta];

            s.coord(s.coord(:,1)== maxC(1) & s.coord(:,3)==minC(3),:) =...
                s.coord(s.coord(:,1)== maxC(1) & s.coord(:,3)==minC(3),:)+[0,0,delta];

            s.coord(s.coord(:,1)== minC(1) & s.coord(:,3)==maxC(3),:) =...
                s.coord(s.coord(:,1)== minC(1) & s.coord(:,3)==maxC(3),:)-[0,0,delta];

            s.coord(s.coord(:,1)== minC(1) & s.coord(:,3)==minC(3),:) =...
                s.coord(s.coord(:,1)== minC(1) & s.coord(:,3)==minC(3),:)+[0,0,delta];

            s.coord(s.coord(:,2)== maxC(2) & s.coord(:,3)==maxC(3),:) =...
                s.coord(s.coord(:,2)== maxC(2) & s.coord(:,3)==maxC(3),:)-[0,0,delta];

            s.coord(s.coord(:,2)== maxC(2) & s.coord(:,3)==minC(3),:) =...
                s.coord(s.coord(:,2)== maxC(2) & s.coord(:,3)==minC(3),:)+[0,0,delta];

            s.coord(s.coord(:,2)== minC(2) & s.coord(:,3)==maxC(3),:) =...
                s.coord(s.coord(:,2)== minC(2) & s.coord(:,3)==maxC(3),:)-[0,0,delta];

            s.coord(s.coord(:,2)== minC(2) & s.coord(:,3)==minC(3),:) =...
                s.coord(s.coord(:,2)== minC(2) & s.coord(:,3)==minC(3),:)+[0,0,delta];

            % COMMENT FOR AIRFOIL
            % s.coord(s.coord(:,1)== maxC(1) & s.coord(:,2)==maxC(2),:) =...
            %     s.coord(s.coord(:,1)== maxC(1) & s.coord(:,2)==maxC(2),:)-[0,delta,0];
            % 
            % s.coord(s.coord(:,1)== maxC(1) & s.coord(:,2)==minC(2),:) =...
            %     s.coord(s.coord(:,1)== maxC(1) & s.coord(:,2)==minC(2),:)+[0,delta,0];
            % 
            % s.coord(s.coord(:,1)== minC(1) & s.coord(:,2)==maxC(2),:) =...
            %     s.coord(s.coord(:,1)== minC(1) & s.coord(:,2)==maxC(2),:)-[0,delta,0];
            % 
            % s.coord(s.coord(:,1)== minC(1) & s.coord(:,2)==minC(2),:) =...
            %     s.coord(s.coord(:,1)== minC(1) & s.coord(:,2)==minC(2),:)+[0,delta,0];

            mS = Mesh.create(s);
        end
        
        function mCoarse = createCoarseMesh(obj)
            s.nsubdomains   = obj.nSubdomains; %nx ny
            s.meshReference = obj.createReferenceCoarseMesh();
            s.tolSameNode   = obj.tolSameNode;
            mRVECoarse      = MeshCreatorFromRVE3D(s);
            [mCoarse,~,~] = mRVECoarse.create();
        end


        function cMesh = createReferenceCoarseMesh(obj)
            mR=obj.referenceMesh;
            xmax = max(mR.coord(:,1));
            xmin = min(mR.coord(:,1));
            ymax = max(mR.coord(:,2));
            ymin = min(mR.coord(:,2));
            zmax = max(mR.coord(:,3));
            zmin = min(mR.coord(:,3));

            % GENERAL CASE
            coord(1,1) = xmin;  coord(1,2) = ymin;   coord(1,3) = zmin;
            coord(2,1) = xmax;  coord(2,2) = ymin;   coord(2,3) = zmin;
            coord(3,1) = xmax;  coord(3,2) = ymax;   coord(3,3) = zmin;
            coord(4,1) = xmin;  coord(4,2) = ymax;   coord(4,3) = zmin;
            coord(5,1) = xmin;  coord(5,2) = ymin;   coord(5,3) = zmax;
            coord(6,1) = xmax;  coord(6,2) = ymin;   coord(6,3) = zmax;
            coord(7,1) = xmax;  coord(7,2) = ymax;   coord(7,3) = zmax;
            coord(8,1) = xmin;  coord(8,2) = ymax;   coord(8,3) = zmax;

            % % Airfoil Abril
            % coord(1,1) = xmin;  coord(1,2) = ymin;   coord(1,3) = zmin;
            % coord(2,1) = xmin;  coord(2,2) = ymin;   coord(2,3) = zmax;
            % coord(3,1) = xmax;  coord(3,2) = ymin;   coord(3,3) = zmax;
            % coord(4,1) = xmax;  coord(4,2) = ymin;   coord(4,3) = zmin;
            % coord(5,1) = xmin;  coord(5,2) = ymax;   coord(5,3) = zmin;
            % coord(6,1) = xmin;  coord(6,2) = ymax;   coord(6,3) = zmax;
            % coord(7,1) = xmax;  coord(7,2) = ymax;   coord(7,3) = zmax;
            % coord(8,1) = xmax;  coord(8,2) = ymax;   coord(8,3) = zmin;

            %Joaquin Airfoil Training
            coord(1,1) = xmin;  coord(1,2) = ymin;   coord(1,3) = zmax;
            coord(2,1) = xmin;  coord(2,2) = ymax;   coord(2,3) = zmax;
            coord(3,1) = xmax;  coord(3,2) = ymax;   coord(3,3) = zmax;
            coord(4,1) = xmax;  coord(4,2) = ymin;   coord(4,3) = zmax;
            coord(5,1) = xmin;  coord(5,2) = ymin;   coord(5,3) = zmin;
            coord(6,1) = xmin;  coord(6,2) = ymax;   coord(6,3) = zmin;
            coord(7,1) = xmax;  coord(7,2) = ymax;   coord(7,3) = zmin;
            coord(8,1) = xmax;  coord(8,2) = ymin;   coord(8,3) = zmin;

            connec = [1 2 3 4 5 6 7 8];    % General case and Joaquin Airfoil
            % connec = [4 1 2 3 8 5 6 7];    % Uncomment for Airfoil Abril
            s.coord = coord;
            s.connec = connec;
            cMesh = Mesh.create(s);  % crea la mesh de 4 nodes
        end

        function createBCapplier(obj)
            s.mesh                  = obj.meshDomain;
            s.boundaryConditions    = obj.boundaryConditions;
            obj.bcApplier           = BCApplier(s);
        end


        function material = createMaterial(obj) 
            switch obj.params.Inclusion
                case {'Hole'}
                    [young,poisson] = obj.computeElasticProperties();
                    s.type          = 'ISOTROPIC';
                    s.ptype         = 'ELASTIC';
                    s.ndim          = obj.meshDomain.ndim;
                    s.young         = young;
                    s.poisson       = poisson;
                    material        = Material.create(s);
                case 'Material'
                    obj.createDesignVariable()
                    m= obj.createMaterialInterpolator();
                    s.type                 = 'DensityBased';
                    s.density              = obj.designVariable;
                    s.materialInterpolator = m;
                    s.dim                  = '3D';
                    s.mesh                 = obj.meshDomain;
                    material = Material.create(s);
                    material.setDesignVariable({obj.designVariable.fun})
                    material = material.obtainTensor();
            end
        end

        function [young,poisson] = computeElasticProperties(obj)
            % GENERAL CASE
            E  = 1;
            nu = 1/3; 

            % AIRFOIL JOAQUIN
            % E  = 70000;
            % nu = 0.3;
            young   = ConstantFunction.create(E,obj.meshDomain);
            poisson = ConstantFunction.create(nu,obj.meshDomain);
        end

        function createDesignVariable(obj)
            ls= obj.computeLevelSet();
            sUm.backgroundMesh = obj.meshDomain;
            sUm.boundaryMesh   = obj.meshDomain.createBoundaryMesh;
            uMesh              = UnfittedMesh(sUm);
            uMesh.compute(-ls);
            % Mprint=UnfittedMesh(sUm);
            % Mprint.compute(ls);
            % obj.unfittedMesh=Mprint;
            funLS        = CharacteristicFunction.create(uMesh);
            s.filterType = 'LUMP';
            s.mesh       = obj.meshDomain;
            s.trial      = LagrangianFunction.create(obj.meshDomain,1,'P1');
            f = Filter.create(s);
            s.fun      = f.compute(funLS,2);
            s.type     = 'Density';
            s.plotting = false;
            dens = DesignVariable.create(s);
            obj.designVariable = dens;
        end

         function ls=computeLevelSet(obj)
            [x0,y0,z0] = obj.computeSubdomainCentroid();
            [Nx,Ny,Nz] = size(obj.r);
            GeomParams(Nx,Ny,Nz) = struct('type',[],'radius',[],'xCoorCenter',[],'yCoorCenter',[],'zCoorCenter',[]);

            for i = 1:Nx
                for j = 1:Ny
                    for k = 1:Nz
                        GeomParams(i,j,k).type        = "Sphere";
                        GeomParams(i,j,k).radius      = obj.r(i,j,k);
                        GeomParams(i,j,k).xCoorCenter = x0(i,j,k);
                        GeomParams(i,j,k).yCoorCenter = y0(i,j,k);
                        GeomParams(i,j,k).zCoorCenter = z0(i,j,k);
                    end
                end
            end

            s.type        = 'GivenPattern';
            s.paramsList  = GeomParams;
            g             = GeometricalFunction(s);
            phiFun        = g.computeLevelSetFunction(obj.meshDomain);
            ls            = phiFun.fValues;
         end

         function [x0,y0,z0]= computeSubdomainCentroid(obj)
            % [Nx, Ny] = size(obj.r);
            xMin = min(obj.meshDomain.coord(:,1));
            xMax = max(obj.meshDomain.coord(:,1));
            yMin = min(obj.meshDomain.coord(:,2));
            yMax = max(obj.meshDomain.coord(:,2));
            zMin = min(obj.meshDomain.coord(:,3));
            zMax = max(obj.meshDomain.coord(:,3));

            dx = obj.xmax-obj.xmin;
            dy = obj.ymax-obj.ymin;
            dz = obj.zmax-obj.zmin;

            x_center = xMin + dx/2 : dx : xMax - dx/2;
            y_center = yMin + dy/2 : dy : yMax - dy/2;
            z_center = zMin + dz/2 : dz : zMax - dz/2;
            [y0, x0, z0] = ndgrid(y_center, x_center, z_center);
            
        end

        function m=createMaterialInterpolator(obj)
            E0 = 1e-3;
            nu0 = 1/3;
            ndim = obj.meshDomain.ndim;
            matA.shear = IsotropicElasticMaterial.computeMuFromYoungAndPoisson(E0,nu0);
            matA.bulk  = IsotropicElasticMaterial.computeKappaFromYoungAndPoisson(E0,nu0,ndim);

            E1 = 1;
            nu1 = 1/3;
            matB.shear = IsotropicElasticMaterial.computeMuFromYoungAndPoisson(E1,nu1);
            matB.bulk  = IsotropicElasticMaterial.computeKappaFromYoungAndPoisson(E1,nu1,ndim);

            s.interpolation  = 'SIMPALL';
            s.dim            = '2D';
            s.matA           = matA;
            s.matB           = matB;

            m = MaterialInterpolator.create(s);
        end

        function [Dir,PL] = createRawBoundaryConditions(obj)
            minx = min(obj.meshDomain.coord(:,1));
            maxx = max(obj.meshDomain.coord(:,1));
            miny = min(obj.meshDomain.coord(:,2));
            maxy = max(obj.meshDomain.coord(:,2));
            minz = min(obj.meshDomain.coord(:,3));
            maxz = max(obj.meshDomain.coord(:,3));
            tolBound = obj.tolSameNode;
            isLeft   = @(coor) (abs(coor(:,1) - minx)   < tolBound);
            isRight  = @(coor) (abs(coor(:,1) - maxx)   < tolBound);
            isBottom = @(coor) (abs(coor(:,3) - minz)   < tolBound);
            
            Dir{1}.domain    = @(coor) isLeft(coor);%| isRight(coor) ;
            Dir{1}.direction = [1,2,3];
            Dir{1}.value     = 0;

            PL.domain    = @(coor) isRight(coor);
            PL.direction = 2;        % 3--> general    2--> Airfoil
            PL.value     = -1;       %Set displacement intensity 
        end 

        function [bc,Dir,PL] = createBoundaryConditions(obj,mesh)
            [Dir,PL]  = obj.createRawBoundaryConditions();
            dirichletFun = [];
            for i = 1:numel(Dir)
                dir = DirichletCondition(obj.meshDomain, Dir{i});
                dirichletFun = [dirichletFun, dir];
            end
            pointload = TractionLoad(mesh,PL,'DIRAC');
            s.pointloadFun = pointload;
            s.dirichletFun = dirichletFun;
            s.periodicFun  =[];
            s.mesh         = mesh;
            bc             = BoundaryConditions(s);
        end


        function [LHSr,RHSr,lhs] = createElasticProblem(obj)
            u = LagrangianFunction.create(obj.meshDomain,obj.meshDomain.ndim,'P1');
            material = obj.createMaterial();
            [lhs,LHSr] = obj.computeStiffnessMatrix(u,obj.meshDomain,material);
            RHSr       = obj.computeForces(lhs,u);
        end

        function [LHS,LHSr] = computeStiffnessMatrix(obj,dispFun,mesh,C)
            f = @(u,v) DDP(SymGrad(v),DDP(C,SymGrad(u)));
            LHS= IntegrateLHS(f,dispFun,dispFun,mesh,'Domain',2);
            LHSr = obj.bcApplier.fullToReducedMatrixDirichlet(LHS);     
        end


        function RHS =  computeForces(obj,stiffness,u)
            bc  = obj.boundaryConditions;
            t   = bc.tractionFun;
            rhs = zeros(u.nDofs,1);
            if ~isempty(t)
                for i = 1:numel(t)
                    rhsi = t(i).computeRHS(u);
                    rhs  = rhs + rhsi;
                end
            end
             bc      = obj.boundaryConditions;
             dirich  = bc.dirichlet_dofs;
             dirichV = bc.dirichlet_vals;
             if ~isempty(dirich)
                 R = -stiffness(:,dirich)*dirichV;
             else
                 R = zeros(sum(u.nDofs(:)),1);
             end
             rhs = rhs+R;
             RHS = obj.bcApplier.fullToReducedVectorDirichlet(rhs);
        end

         function Meifem = createEIFEMPreconditioner(obj,dir,iC,lG,bS,iCR,dMesh)
            mR = obj.referenceMesh;
            EIFEMfilename = obj.fileNameEIFEM;
            s.RVE           = TrainedRVE(EIFEMfilename);
            s.mesh          = obj.createCoarseMesh();
            s.DirCond       = dir;
            s.nSubdomains   = obj.nSubdomains;
            eifem           = EIFEM(s);
            
            ss.ddDofManager = obj.createDomainDecompositionDofManager(iC,lG,bS,mR,iCR);
            ss.EIFEMsolver  = eifem;
            ss.bcApplier    = obj.bcApplier;
            ss.dMesh        = dMesh;
            ss.type         = 'EIFEM';
            eP = Preconditioner.create(ss);
            Meifem = @(r) eP.apply(r);
        end        


        function Mcoarse = createCoarseNNPreconditioner(obj,mR,dir,iC,lG,bS,iCR,dMesh)
            p=obj.params;
            meshName    =  p.nelem+"x"+p.nelem;
            switch p.Option
                case 'Dataset'
                    nameFile=obj.computeNameFile();
                    obj.loadDataset(nameFile);
                case 'NN'
                    filePath = fullfile("AbrilTFGfiles","Data/Sphere/",p.Training,"K_NN.mat");
                    load(filePath,"K_NN");
                    filePath = fullfile("AbrilTFGfiles","Data/Sphere/",p.Training,"T_NN.mat");
                    load(filePath,"T_NN","pol_deg");
            end

            RVE = cell(obj.nSubdomains(1,2),obj.nSubdomains(1,1));

            for i = 1:obj.nSubdomains(1,2)
                for j = 1:obj.nSubdomains(1,1)
                    RVE{i,j}.ndimf = 2;

                    switch p.Option
                        case 'Dataset'
                            RVE{i,j}.Kcoarse= obj.data.K{i,j};
                            RVE{i,j}.U= obj.data.T{i,j}; 
                        case 'NN'
                            RVE{i,j}.Kcoarse = computeKcoarse_NN(K_NN,obj.r(i,j));
                            RVE{i,j}.U       = computeT_NN(obj.cellMesh{i,j},obj.r(i,j),T_NN,pol_deg);
                        case 'Hybrid'
                            RVE{i,j}.Kcoarse = computeKcoarse_NN(K_NN,obj.r(i,j));
                            RVE{i,j}.U       = computeT_Hybrid(basis,obj.r(i,j),Q_NN,pol_deg);
                    end
                end
            end

            s.RVE           = RVE;
            s.mesh          = obj.createCoarseMesh();
            s.DirCond       = dir;
            s.nSubdomains   = obj.nSubdomains;
            coarseSolver    = Coarse(s);

            ss.ddDofManager = obj.createDomainDecompositionDofManager(iC,lG,bS,mR,iCR);
            ss.Coarsesolver = coarseSolver;
            ss.bcApplier = obj.bcApplier;
            ss.dMesh     = dMesh;
            ss.type = 'Coarse';
            eP = Preconditioner.create(ss);
            Mcoarse = @(r) eP.apply(r);

            obj.ddDofManager=ss.ddDofManager;
            close all % Ho he afegit pq s'obren fig sense info
        end

        function NameFile=computeNameFile(obj)
            n=obj.params.nelem;
            rad=obj.r;
            meshName=n+"x"+n;
            name=strings(size(rad,1),size(rad,2));
            for i=1:size(rad,1)
                for j=1:size(rad,2)
                    name(i,j) = strrep("r"+num2str(rad(i,j), '%.4f'), ".", "_")+"-"+meshName+".mat";
                end
            end
            NameFile=name;
        end


        function loadDataset(obj,name)
            p=obj.params;
            n=p.nelem;
            Taux=cell(size(name,1),size(name,2));
            Kaux=cell(size(name,1),size(name,2));
            meshName=n+"x"+n;
            for i=1:size(name,1)
                for j=1:size(name,2)
                    switch p.Inclusion
                        case {'Material','HoleRaul'}
                                filePath = fullfile("AbrilTFGfiles","Data/Sphere/",p.Training,meshName,name(i,j));
                        case 'Hole'
                                filePath = fullfile('AbrilTFGfiles', 'Data/Sphere/',p.Training,'hole',name(i,j));
                    end
                    load(filePath,"T","Kcoarse");
                    Taux{i,j}=T;
                    Kaux{i,j}=Kcoarse;
                end
            end
            obj.data.T=Taux;
            obj.data.K=Kaux;
        end


        function d = createDomainDecompositionDofManager(obj,iC,lG,bS,mR,iCR)
            s.nSubdomains     = obj.nSubdomains;
            s.interfaceConnec = iC;
            s.interfaceConnecReshaped = iCR;
            s.locGlobConnec   = lG;
            s.nBoundaryNodes  = bS{1}.mesh.nnodes;
            s.nReferenceNodes = mR.nnodes;
            s.nNodes          = obj.meshDomain.nnodes;
            s.nDimf           = obj.meshDomain.ndim;
            d = DomainDecompositionDofManager3D(s);
        end

        function Milu = createILUpreconditioner(~,LHS)
            s.LHS = sparse(LHS);
            s.type = 'ILU';
            M = Preconditioner.create(s);
            Milu = @(r) M.apply(r);
        end

    end

    methods (Static, Access = public)

        function J = computeTotalEnergy(x,A,b)
            J = 0.5*x'*A(x)-b'*x;
        end

    end

end
