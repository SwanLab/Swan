classdef MainSteadyConvectionDiffusion1D < handle

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
            obj.a  = 1;
            obj.nu = 0.01;
        end

        function createMesh(obj)
            nEl           = 10;
            xnode         = 0:1/nEl:1; % Podríem demanar com a 1a pregunta que modifiquin els inputs amb Peclet
            s.coord       = xnode';
            s.connec(:,1) = 1:length(xnode)-1;
            s.connec(:,2) = 2:length(xnode);
            obj.mesh      = Mesh.create(s);
        end

        function createFiniteElementFunction(obj)
            order     = 'P1'; % P1/P2
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
            s.mesh       = obj.mesh;
            s.trial      = obj.trial;
            s.sHandle    = obj.source;
            s.dirValues  = obj.dirValues;
            s.stab       = 'Galerkin';
            prob         = SteadyConvectionDiffusionProblem(s);
            obj.solution = prob.compute(obj.a,obj.nu);
        end
    end
end
