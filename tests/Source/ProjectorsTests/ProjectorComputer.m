classdef ProjectorComputer < handle

    properties (Access = public)
        computation
    end

    properties (Access = private)
        problemType
        testName
        problemTestName
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
        end

        function computeProjection(obj)
            switch obj.problemType
                case {'FEM'}
                    strain = squeeze(obj.computation.variables.strain)';
                    z.fValues = strain;
                    strainFun = P0Function(z);
                    a.fileName = obj.problemTestName;
                    s = FemDataContainer(a);
                    bb.mesh   = s.mesh;
                    bb.connec = s.mesh.connec;
                    projector = Projector_P0toP1(bb);
                    p1strain = projector.project(strainFun);
                    nrows = numel(p1strain.fValues);
                    obj.computation.variables.xP = reshape(p1strain.fValues,[nrows,1]);
                case {'TOPOPT'}
                    %...
            end
        end

    end
end