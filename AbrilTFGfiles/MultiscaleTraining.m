classdef MultiscaleTraining < handle
       
    properties (Access = private)
        mesh
        material
    end

    properties  (Access = private)
        boundaryMeshJoined
        localGlobalConnecBd
    end


    methods (Access = public)

        function [obj, u, L, Kc] = MultiscaleTraining(cParams)
            obj.init(cParams)
            [obj.boundaryMeshJoined, obj.localGlobalConnecBd] = obj.mesh.createSingleBoundaryMesh();           
            [u,L,Kc] = obj.solveElasticHarmonicExtenstion();
        end
    end

    methods (Access = private)

        function init(obj,cParams)
            obj.mesh     = cParams.mesh;
            obj.material = cParams.material;
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
            s.uFun         = LagrangianFunction.create(obj.mesh, obj.mesh.ndim, 'P1');
            s.lambdaFun    = LagrangianFunction.create(obj.boundaryMeshJoined, obj.mesh.ndim, 'P1'); 
            s.material      = obj.material;   
            s.dirichletFun = obj.createDirichletFunctions();
            s.localGlobalConnecBd = obj.localGlobalConnecBd;
            e  = ElasticHarmonicExtension(s);
            [u,L,K] = e.solve();
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
