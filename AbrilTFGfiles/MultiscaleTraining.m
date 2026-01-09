classdef MultiscaleTraining < handle
    
    properties (Access = public)
        bcApplier
        centroids
        nelem
    end
    
    properties (Access = private)
        stiffness
        radius
        nodeDirection
        physicalProblem
        boundaryConditions
        problemSolver
        forces
        uFun
        strainFun 
        Inclusion
    end

    properties  (Access = protected)
        mesh
        material
        displacementFun
        boundaryMesh
        boundaryMeshJoined
        localGlobalConnecBd

    end


    methods (Access = public)

        function [obj, K, u, L, mesh,Kcoarse] = MultiscaleTraining(meshRef,r,params)
            obj.init(meshRef,r,params)

            obj.boundaryMesh = obj.mesh.createBoundaryMesh();
            [obj.boundaryMeshJoined, obj.localGlobalConnecBd] = obj.mesh.createSingleBoundaryMesh();

            
            [u, L,K] = obj.solveElasticHarmonicExtenstion();
            mesh = obj.mesh;
            
            %K=obj.stiffness;
            Kcoarse=u.'*K*u;
            
        end
    end


    methods (Access = private)

        function init(obj,meshRef,r,params)
            % close all;
            % clc;
            obj.radius = r;
            obj.nodeDirection = 1;
            obj.nelem=params.nelem;
            obj.Inclusion=params.Inclusion;
            obj.mesh=meshRef;
        end

      %  function createMesh(obj)
      %       obj.boundaryMesh = obj.mesh.createBoundaryMesh();
      %      [obj.boundaryMeshJoined, obj.localGlobalConnecBd] = obj.mesh.createSingleBoundaryMesh();
      %  end

       % function mesh = createReferenceMesh(obj)
       %     n       = obj.nelem;
       %     x1      = linspace(-1,1,n);
       %     x2      = linspace(-1,1,n);
       %     [xv,yv] = meshgrid(x1,x2);
       %     [F,V]   = mesh2tri(xv,yv,zeros(size(xv)),'x');
       %     s.coord  = V(:,1:2);
       %     s.connec = F;
       %     mesh = Mesh.create(s);
       % end

        function levelSet = createLevelSetFunction(obj,bgMesh)
            sLS.type        = 'CircleInclusion';
            sLS.xCoorCenter = 0;
            sLS.yCoorCenter = 0;
            sLS.radius      = obj.radius;
            g               = GeometricalFunction(sLS);
            lsFun           = g.computeLevelSetFunction(bgMesh);
            levelSet        = lsFun.fValues;
        end

        function uMesh = computeUnfittedMesh(~,bgMesh,levelSet)
            sUm.backgroundMesh = bgMesh;
            sUm.boundaryMesh   = bgMesh.createBoundaryMesh();
            uMesh              = UnfittedMesh(sUm);
            uMesh.compute(levelSet);
        end

        function [u,L,K] = solveElasticHarmonicExtenstion(obj)
            s.mesh         = obj.mesh;
            s.boundaryMesh = obj.boundaryMeshJoined;
            s.uFun      = LagrangianFunction.create(obj.mesh, obj.mesh.ndim, 'P1');
            s.lambdaFun = LagrangianFunction.create(obj.boundaryMeshJoined, obj.mesh.ndim, 'P1'); 
            s.material  = obj.createMaterial();   
            s.dirichletFun = obj.createDirichletFunctions();
            s.localGlobalConnecBd = obj.localGlobalConnecBd;
            e  = ElasticHarmonicExtension(s);
            [u,L,K] = e.solve();
        end


        function material = createMaterial(obj)
            [young,poisson] = obj.computeElasticProperties(obj.mesh);
            s.type       = 'ISOTROPIC';
            s.ptype      = 'ELASTIC';
            s.ndim       = obj.mesh.ndim;
            s.young      = young;
            s.poisson    = poisson;
            material = Material.create(s);            
        end

        function [young,poisson] = computeElasticProperties(obj,mesh)
            E  = 1;
            nu = 1/3;
            r  = obj.radius;
            switch obj.Inclusion
                case {'Hole','HoleRaul'}
                    young   = ConstantFunction.create(E,mesh);
                    poisson = ConstantFunction.create(nu,mesh);
                case 'Material'
                    E2 = E/1000;
                    x0=mean(mesh.coord(:,1));
                    y0=mean(mesh.coord(:,2));
                    f   = @(x) (sqrt((x(1,:,:)-x0).^2+(x(2,:,:)-y0).^2)<r)*E2 + ...
                                (sqrt((x(1,:,:)-x0).^2+(x(2,:,:)-y0).^2)>=r)*obj.E ; 
                    young   = AnalyticalFunction.create(f,mesh);
                    poisson = ConstantFunction.create(nu,mesh);
            end          
        end

        function createBCApplyerHere(obj, cParams)
            s.mesh = obj.mesh;
            s.boundaryConditions = cParams.boundaryConditions;
            obj.bcApplier = BCApplier(s);
        end

        function uD = createDirichletFunctions(obj)
            f1x = @(x) [1/(4)*(1-x(1,:,:)).*(1-x(2,:,:));...
                    0*x(2,:,:)  ];
            f2x = @(x) [1/(4)*(1+x(1,:,:)).*(1-x(2,:,:));...
                    0*x(2,:,:)  ];
            f3x = @(x) [1/(4)*(1+x(1,:,:)).*(1+x(2,:,:));...
                    0*x(2,:,:)  ];
            f4x = @(x) [1/(4)*(1-x(1,:,:)).*(1+x(2,:,:));...
                    0*x(2,:,:)  ];

            f1y = @(x) [0*x(1,:,:);...
                    1/(4)*(1-x(1,:,:)).*(1-x(2,:,:))  ];
            f2y = @(x) [0*x(1,:,:);...
                    1/(4)*(1+x(1,:,:)).*(1-x(2,:,:))  ];
            f3y = @(x) [0*x(1,:,:);...
                    1/(4)*(1+x(1,:,:)).*(1+x(2,:,:))  ];
            f4y = @(x) [0*x(1,:,:);...
                    1/(4)*(1-x(1,:,:)).*(1+x(2,:,:))  ];

            fB     = {f1x f1y f2x f2y f3x f3y f4x f4y}; 
            nfun = size(fB,2);
            uD=cell(1,8);
            for i=1:nfun
                uD{i}  = AnalyticalFunction.create(fB{i},obj.boundaryMeshJoined);
            end
        end

    end
    
end
