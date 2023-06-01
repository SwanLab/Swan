classdef ProjectorComputer < handle

    properties (Access = public)
        computation
        variables
    end

    properties (Access = private)
        problemTestName
        projectTo
        fun
        funProj
        mesh
        var
    end

    properties (Access = private)
        problemType
        testName
    end

    methods (Access = public)

        function obj = ProjectorComputer(cParams)
            obj.init(cParams);
            obj.saveConfiguration();
        end

        function compute(obj)
            s.testName = obj.problemTestName;
            computer   = TestComputer.create(obj.problemType, s);
            computer.compute();
            obj.computation = computer.computation;
            obj.computeProjection();
            nrows = numel(obj.funProj.fValues);
            obj.variables.xP = reshape(obj.funProj.fValues,[nrows,1]);
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.problemType = cParams.problemType;
            obj.testName    = cParams.testName;
        end

        function saveConfiguration(obj)
            run(obj.testName);
            obj.problemTestName = problemToSolve;
            obj.projectTo       = projectDestination;
            obj.var             = variable;
        end

        function computeProjection(obj)
            switch obj.problemType
                case {'FEM'}
                    obj.readMesh();
                    obj.selectOrigin();
                    obj.selectProjector();
                case {'TOPOPT'}
                    %...
            end
        end

        function readMesh(obj)
            a.fileName = obj.problemTestName;
            s          = FemDataContainer(a);
            obj.mesh   = s.mesh;
        end

        function selectOrigin(obj)
            obj.fun = obj.computation.(obj.var);
        end

        function selectProjector(obj)
            s.mesh          = obj.mesh;
            s.connec        = obj.mesh.connec;
            s.projectorType = obj.projectTo;
            projector       = Projector.create(s);
            obj.funProj     = projector.project(obj.fun);
        end

    end
end