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
            obj.variables.xP = obj.funProj.fValues;
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
                    obj.project();
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

        function project(obj)
            obj.funProj = obj.fun.project(obj.projectTo);
        end

    end
end