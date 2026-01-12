classdef MultiscaleTraining < handle
       
    properties (Access = private)
        mesh
        material
    end

    methods (Access = public)

        function obj = MultiscaleTraining(cParams)
            obj.init(cParams);
        end

        function [u,L,Kc] = train(obj)
            [bMesh,lGCBd]   = obj.mesh.createSingleBoundaryMesh();
            s.mesh          = obj.mesh;
            s.uFun          = LagrangianFunction.create(obj.mesh, obj.mesh.ndim, 'P1');
            s.lambdaFun     = LagrangianFunction.create(bMesh,obj.mesh.ndim, 'P1');
            s.material      = obj.material;
            s.dirichletFun  = obj.createDirichletFunctions(bMesh);
            s.localGlobalConnecBd = lGCBd;
            e  = ElasticHarmonicExtension(s);
            [u,L,Kc] = e.solve();
        end
    end

    methods (Access = private)

        function init(obj,cParams)
            obj.mesh     = cParams.mesh;
            obj.material = cParams.material;
        end

        function uD = createDirichletFunctions(obj,bMesh)
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
                uD{i}  = AnalyticalFunction.create(fB{i},bMesh);
            end
        end

    end
    
end
