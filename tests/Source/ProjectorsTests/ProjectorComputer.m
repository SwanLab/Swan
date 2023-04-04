classdef ProjectorComputer < handle

    properties (Access = public)
        computation
    end

    properties (Access = private)
        problemTestName
        projectFrom
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
            obj.projectFrom     = origin;
            obj.projectTo       = projectDestination;
            obj.var             = variable;
        end

        function computeProjection(obj)
            switch obj.problemType
                case {'FEM'}
                    obj.readMesh();
                    obj.selectOrigin();
                    obj.selectProjector();
                    nrows = numel(obj.funProj.fValues);
                    obj.computation.variables.xP = reshape(obj.funProj.fValues,[nrows,1]);
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
            val = obj.computation.variables.(obj.var);
            switch obj.projectFrom
                case {'P0'}
                    val = squeeze(val)';
                case {'P1'}
                    val = reshape(val,[obj.mesh.ndim,obj.mesh.nnodes])';
            end
            z.fValues      = val;
%             z.connec       = obj.mesh.connec;
%             z.type         = obj.mesh.type;
            z.mesh         = obj.mesh;
            z.functionType = obj.projectFrom;
            obj.fun        = FeFunction.create(z);
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