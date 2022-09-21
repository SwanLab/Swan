classdef ProjectorComputer < handle

    properties (Access = public)
        computation
    end

    properties (Access = private)
        problemType
        testName
        problemTestName
        projectFrom
        projectTo
        fun
        funProj
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
        end

        function computeProjection(obj)
            switch obj.problemType
                case {'FEM'}
                    a.fileName = obj.problemTestName;
                    s = FemDataContainer(a);
                    obj.selectOrigin(s);
                    obj.selectProjector(s);
                    nrows = numel(obj.funProj.fValues);
                    obj.computation.variables.xP = reshape(obj.funProj.fValues,[nrows,1]);
                case {'TOPOPT'}
                    %...
            end
        end

        function selectOrigin(obj,s)
            switch obj.projectFrom
                case {'P0'}
                    strain = squeeze(obj.computation.variables.strain)';
                    z.fValues = strain;
                    obj.fun = P0Function(z);
                case {'P1'}
                    uCol = obj.computation.variables.d_u;
                    u    = reshape(uCol,[s.mesh.ndim,s.mesh.nnodes])';
                    z.connec = s.mesh.connec;
                    z.type   = s.mesh.type;
                    z.fNodes = u;
                    obj.fun  = P1Function(z);
            end
        end

        function selectProjector(obj,s)
            switch obj.projectTo
                case {'P0'}
                    b.mesh   = s.mesh;
                    b.connec = s.mesh.connec;
                    b.nelem  = size(s.mesh.connec,1);
                    b.nnode  = size(s.mesh.connec,2);
                    b.npnod  = size(s.mesh.coord,1);
                    projector = Projector_P1toP0(b);
                    obj.funProj = projector.project(obj.fun);
                case {'P1'}
                    b.mesh   = s.mesh;
                    b.connec = s.mesh.connec;
                    projector = Projector_P0toP1(b);
                    obj.funProj = projector.project(obj.fun);
            end
        end

    end
end