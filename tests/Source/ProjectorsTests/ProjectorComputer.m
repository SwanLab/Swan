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
            switch obj.projectFrom
                case {'P0'}
                    varP0 = squeeze(obj.computation.variables.strain)';
                    z.fValues = varP0;
                    z.connec  = obj.mesh.connec;
                    z.type    = obj.mesh.type;
                    obj.fun   = P0Function(z);
                case {'P1'}
                    varP1 = obj.computation.variables.d_u;
                    varP1 = reshape(varP1,[obj.mesh.ndim,obj.mesh.nnodes])';
                    z.connec  = obj.mesh.connec;
                    z.type    = obj.mesh.type;
                    z.fValues = varP1;
                    obj.fun   = P1Function(z);
            end
        end

        function selectProjector(obj)
            switch obj.projectTo
                case {'P0'}
                    b.mesh      = obj.mesh;
                    b.connec    = obj.mesh.connec;
                    b.nelem     = size(obj.mesh.connec,1);
                    b.nnode     = size(obj.mesh.connec,2);
                    b.npnod     = size(obj.mesh.coord,1);
                    b.projectorType = 'toP0';
                    projector   = Projector.create(b);
                    obj.funProj = projector.project(obj.fun);
                case {'P1'}
                    b.mesh      = obj.mesh;
                    b.connec    = obj.mesh.connec;
                    b.projectorType = 'toP1';
                    projector   = Projector.create(b);
                    obj.funProj = projector.project(obj.fun);
                case {'P1disc'}
                    b.mesh      = obj.mesh;
                    b.connec    = obj.mesh.connec;
                    projector   = Projector_toP1Discontinuous(b);
                    obj.funProj = projector.project(obj.fun);
            end
        end

    end
end