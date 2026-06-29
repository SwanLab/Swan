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
        tFrame
        tCross
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
            % Usol   = LHS\RHS;
            % Ufull  = obj.bcApplier.reducedToFullVectorDirichlet(Usol); 
            t_direct=toc

            % PRECONDITIONERS
            Milu        = obj.createILUpreconditioner(LHS);
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
            % 
            % tic %SOLVE THE CASE WITH STANDARD CG
            % [~,obj.residualCG,errCG, errCG] = PCG.solve(LHSf,RHSf,x0,Mid,tol,Usol,obj.meshDomain,obj.bcApplier);
            % t_CG=toc
            % tic  % SOLVE THE CASE WITH CG+ ILU
            % [~,obj.residualILU,errILU, errAnormILU] = PCG.solve(LHSf,RHSf,x0,Milu,tol,Usol,obj.meshDomain,obj.bcApplier);
            % t_ILU=toc
            tic %SOLVE THE CASE WITH PRECONDITIONING ILU+EIFEM+ILU
            % [uPCG,obj.residualPCG,obj.errPCG,obj.errAnormPCG] = PCG.solve(LHSf,RHSf,x0,Mmult,tol,Usol,obj.meshDomain,obj.bcApplier);
            [uPCG,obj.residualPCG,obj.errPCG,obj.errAnormPCG] = PCG.solve(LHSf,RHSf,x0,Mmult,tol);
            t_PCG=toc
            xFull = obj.bcApplier.reducedToFullVectorDirichlet(uPCG);
            

            uDomain = obj.bcApplier.reducedToFullVectorDirichlet(uPCG);
            % uDomain = obj.ddDofManager.global2local(uDomain);
            
            % LAGRANGIAN FUN SOLUTIONS
            % s.mesh     = obj.meshDomain;
            % s.ndimf    = obj.meshDomain.ndim;
            % s.order    = 'P1';
            % s.fValues  = reshape(xFull,3,[])';
            % obj.Sol    = LagrangianFunction(s); %Preconditioned sol  
            % s.fValues = reshape(Ufull,3,[])';
            % obj.SolExact=LagrangianFunction(s); %Exact sol
            % 
            % obj.print(xFull,"Prova1stIterInnerUniformCube");
            % obj.print(Ufull,"ProvaSolExactInnerUniformCube");
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
            p.Training      = cParams.Training;    % 'EIFEM'/'Multiscale'/'EIFisol'
            p.Inclusion     = cParams.Inclusion;   % 'Hole'/'Material'/'HoleRaul'   --> Hole: just for constant r
            p.Option        = cParams.Option;      % 'Dataset'/'NN'/'Direct'
            p.nelem         = cParams.nelem;       %  Mesh refining
            p.Geometry      = cParams.Geometry;    % 'Sphere'/'Cube'
            obj.params      = p;
            if isfield(cParams,'r')
                obj.r       = cParams.r;
                [Ny,Nx,Nz] = size(obj.r);
            elseif isfield(cParams,'tFrame')
                obj.tFrame  = cParams.tFrame;
                obj.tCross  = cParams.tCross;
                [Ny,Nx,Nz] = size(obj.tFrame);
            end
            obj.nSubdomains = [Nx Ny Nz];
            % obj.nSubdomains = [35 1 1];    % UNCOMMENT JUST FOR AIRFOIL
            obj.tolSameNode = 1e-11;   % 1E-10--> general case   1E-6 --> airfoil
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

            delta=1E-9;

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
            % coord(1,1) = xmin;  coord(1,2) = ymin;   coord(1,3) = zmax;
            % coord(2,1) = xmin;  coord(2,2) = ymax;   coord(2,3) = zmax;
            % coord(3,1) = xmax;  coord(3,2) = ymax;   coord(3,3) = zmax;
            % coord(4,1) = xmax;  coord(4,2) = ymin;   coord(4,3) = zmax;
            % coord(5,1) = xmin;  coord(5,2) = ymin;   coord(5,3) = zmin;
            % coord(6,1) = xmin;  coord(6,2) = ymax;   coord(6,3) = zmin;
            % coord(7,1) = xmax;  coord(7,2) = ymax;   coord(7,3) = zmin;
            % coord(8,1) = xmax;  coord(8,2) = ymin;   coord(8,3) = zmin;

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
            Mprint=UnfittedMesh(sUm);
            Mprint.compute(ls);
            aa=Mprint.createInnerMesh;
            obj.unfittedMesh=Mprint;
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
            Nx=obj.nSubdomains(2);
            Ny=obj.nSubdomains(1);
            Nz=obj.nSubdomains(3);
            GeomParams(Nx,Ny,Nz) = struct('type',[],'radius',[],'xCoorCenter',[],'yCoorCenter',[],'zCoorCenter',[],'tFrame',[],'tCross',[]);

            switch obj.params.Geometry
                case 'Sphere'
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

                case 'Cube'
                    for i = 1:Nx
                        for j = 1:Ny
                            for k = 1:Nz
                                GeomParams(i,j,k).type        = "Cube";
                                GeomParams(i,j,k).length      = obj.r(i,j,k);
                                GeomParams(i,j,k).xCoorCenter = x0(i,j,k);
                                GeomParams(i,j,k).yCoorCenter = y0(i,j,k);
                                GeomParams(i,j,k).zCoorCenter = z0(i,j,k);
                            end
                        end
                    end

                case 'Lattice3D'
                    for i = 1:Nx
                        for j = 1:Ny
                            for k = 1:Nz
                                GeomParams(i,j,k).type        = "CrossedSquare3D";
                                GeomParams(i,j,k).length      = 2;
                                GeomParams(i,j,k).tFrame      = obj.tFrame(i,j,k);
                                GeomParams(i,j,k).tCross      = obj.tCross(i,j,k);
                                GeomParams(i,j,k).xCoorCenter = x0(i,j,k);
                                GeomParams(i,j,k).yCoorCenter = y0(i,j,k);
                                GeomParams(i,j,k).zCoorCenter = z0(i,j,k);
                            end
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
            s.dim            = '3D';
            s.matA           = matA;
            s.matB           = matB;

            m = MaterialInterpolator.create(s);
        end

        function [Dir,PL] = createRawBoundaryConditions(obj)
            xMin    = min(obj.meshDomain.coord(:,1));
            xMax    = max(obj.meshDomain.coord(:,1));
            yMax    = max(obj.meshDomain.coord(:,2));
            yMin    = min(obj.meshDomain.coord(:,2));
            zMin    = min(obj.meshDomain.coord(:,3));
            zMax    = max(obj.meshDomain.coord(:,3));
            Ly      = yMax-yMin;
            Lx      = xMax-xMin;
            Lz      = zMax-zMin;
            tolBound = obj.tolSameNode;
            isLeft   = @(coor) (abs(coor(:,1) - xMin)   < tolBound);
            isRight  = @(coor) (abs(coor(:,1) - zMax)   < tolBound);
            isBottom = @(coor) (abs(coor(:,3) - zMin)   < tolBound);

            Dir{1}.domain    = @(coor) isLeft(coor);%| isRight(coor) ;
            Dir{1}.direction = [1,2,3];
            Dir{1}.value     = 0;

            PL.domain    = @(coor) isRight(coor);
            PL.direction = 3;        % 3--> general    2--> Airfoil
            PL.value     = -1;       %Set displacement intensity 

            isDir   = @(coor)  abs(coor(:,1)-xMin) < 1e-9;
            isForceDirac = @(coor)  (abs(coor(:,1)-xMax))< 2e-9 &...
                (coor(:,2)>=yMin+0.46* Ly & coor(:,2)<=yMin+0.54* Ly) & ...
                (coor(:,3)>=zMin+0.46* Lz & coor(:,3)<=zMin+0.54* Lz);
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
                    filePath = fullfile("AbrilTFGfiles","Data",p.Geometry,p.Training,"K_NN.mat");
                    load(filePath,"K_NN");
                    filePath = fullfile("AbrilTFGfiles","Data",p.Geometry,p.Training,"T_NN.mat");
                    load(filePath,"T_NN","pol_deg");
            end
            Nx=obj.nSubdomains(2);
            Ny=obj.nSubdomains(1);
            Nz=obj.nSubdomains(3);
            RVE = cell(Nx,Ny,Nz);

            % U=obj.computeT_NN2(T_NN,pol_deg);

            for i = 1:Nx
                for j = 1:Ny
                    for k = 1:Nz
                        RVE{i,j,k}.ndimf = 3;

                        switch p.Option
                            case 'Dataset'
                                RVE{i,j,k}.Kcoarse= obj.data.K{i,j,k};
                                RVE{i,j,k}.U= obj.data.T{i,j,k};
                            case 'NN'
                                if ~isempty(obj.r)
                                    RVE{i,j,k}.Kcoarse = computeKcoarse_NN(K_NN,obj.r(i,j,k),24);
                                    RVE{i,j,k}.U       = computeT_NN(obj.referenceMesh,obj.r(i,j,k),T_NN,pol_deg);
                                elseif ~isempty(obj.tFrame)
                                    RVE{i,j,k}.Kcoarse = computeKcoarse_NN(K_NN,[obj.tFrame(i,j,k),obj.tCross(i,j,k)],24);
                                    RVE{i,j,k}.U       = computeT_NN(obj.referenceMesh,[obj.tFrame(i,j,k),obj.tCross(i,j,k)],T_NN,pol_deg);
                                end
                        end
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
            meshName=n+"x"+n;
            Nx=obj.nSubdomains(2);
            Ny=obj.nSubdomains(1);
            Nz=obj.nSubdomains(3);
            name=strings(Nx,Ny,Nz);
            
            if ~isempty(obj.r)                
                for i = 1:Nx
                    for j = 1:Ny
                        for k = 1:Nz
                            name(i,j,k) = strrep("r"+num2str(obj.r(i,j,k), '%.4f'), ".", "_")+"-"+meshName+".mat";
                        end
                    end
                end
            elseif ~isempty(obj.tFrame)
                t1=obj.tFrame;
                t2=obj.tCross;
                for i = 1:Nx
                    for j = 1:Ny
                        for k = 1:Nz
                            name(i,j,k) = strrep("t1_"+num2str(t1(i,j,k), '%.2f'), ".", "_")+strrep("_t2_"+num2str(t2(i,j,k), '%.2f'), ".", "_")+"-"+meshName+".mat";
                        end
                    end
                end
            end
            NameFile=name;
        end


        function loadDataset(obj,name)
            p=obj.params;
            n=p.nelem;
            [Nx,Ny,Nz] = size(name);
            Taux=cell(Nx,Ny,Nz);
            Kaux=cell(Nx,Ny,Nz);
            meshName=n+"x"+n;

            for i=1:Nx
                for j=1:Ny
                    for k=1:Nz
                        filePath = fullfile("AbrilTFGfiles","Data",p.Geometry,p.Training,meshName,name(i,j,k));
                        load(filePath,"T","Kcoarse");
                        Taux{i,j,k}=T;
                        Kaux{i,j,k}=Kcoarse;
                    end
                end
            end
            obj.data.T=Taux;
            obj.data.K=Kaux;
        end

        function T_trained=computeT_NN2(obj,T_NN,pol_deg)
            coords= obj.referenceMesh.coord;
            sz=size(obj.r);
            N= numel(obj.r);
            nnodes=obj.referenceMesh.nnodes;
            ndim =obj.referenceMesh.ndim;

            r_all= obj.r(:);
            r_in    = repelem(r_all,nnodes,1);
            coords_in = repmat(coords,N,1);

            dataInput= [r_in,coords_in];
            dataFull = Data.buildModel(dataInput,pol_deg);

            T_all  = T_NN.computeOutputValues(dataFull); % 
            T_aux1 = reshape(T_all, ndim, [], nnodes, N);
    
            nModes = size(T_aux1,2);

            T_cells = cell(sz);

            for idx = 1:N

                T_sub = T_aux1(:,:,:,idx);   % ndim × nModes × nNodes

                % colocar como (x1;y1;z1;x2;...)
                T_sub = permute(T_sub, [1 3 2]);   % ndim × nNodes × nModes

                T_cells{idx} = reshape(T_sub, [], nModes);
            end

        end

        function K=computeKcoarse_NN2(obj,K_NN,r)
            K_aux1=K_NN.computeOutputValues(r);
            K_aux2=zeros(8);
            idx=1;
            for n=1:8
                for m=n:8
                    K_aux2(n,m)=K_aux1(idx);
                    idx=idx+1;
                end
            end
            K=K_aux2+triu(K_aux2,1).';
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
