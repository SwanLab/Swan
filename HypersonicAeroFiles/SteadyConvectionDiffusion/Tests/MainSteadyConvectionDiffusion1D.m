classdef MainSteadyConvectionDiffusion1D < handle % Transform to mlx

    properties (Access = private)
        a
        nu
        mesh
        trial
        dirValues
        source
        solution
    end

    methods (Access = public)

        function obj = MainSteadyConvectionDiffusion1D()
            obj.init();
            obj.createMesh();
            obj.createFiniteElementFunction();
            obj.createDirichletConditions();
            obj.createSourceFunction();
            obj.solveProblem();
        end

    end

    methods (Access = private)

        function init(obj)
            obj.a  = 10;
            obj.nu = 0.01;
        end

        function createMesh(obj)
            nEl           = 100;
            xnode         = 0:1/nEl:1;
            s.coord       = xnode';
            s.connec(:,1) = 1:length(xnode)-1;
            s.connec(:,2) = 2:length(xnode);
            obj.mesh      = Mesh.create(s);
        end

        function createFiniteElementFunction(obj)
            order     = 'P1';
            obj.trial = LagrangianFunction.create(obj.mesh,1,order);
        end

        function createDirichletConditions(obj)
            dir.left        = 0;
            dir.right       = 1;
            obj.dirValues   = dir;
        end

        function createSourceFunction(obj)
            obj.source = @(x) 10*exp(-5*x(1,:,:))-4*exp(-x(1,:,:));
        end

        function solveProblem(obj)
            cases = {'Galerkin','Upwind','SUPG'};
            for i=1:length(cases)
                s.mesh       = obj.mesh;
                s.trial      = obj.trial;
                s.sHandle    = obj.source;
                s.dirValues  = obj.dirValues;
                s.stab       = cases{i};
                prob         = SteadyConvectionDiffusionProblem(s);
                obj.solution = prob.compute(obj.a,obj.nu);
                % prob.plot
                hold on
            end
        end
    end
end
