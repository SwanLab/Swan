classdef RHSintegrator_Composite < handle

    properties (GetAccess = public, SetAccess = private)
        integrators
        nInt
        dofs
    end

    properties (Access = private)
        RHScells
        RHSsubcells
        unfittedMesh
        testClass
        test
    end

    methods (Access = public)

        function obj = RHSintegrator_Composite(cParams)
            obj.init(cParams);
            obj.createIntegrators(cParams);
        end

        function f = integrate(obj,nodalFunc)
            f = cell(1,obj.nInt);
            for iInt = 1:obj.nInt
                f{iInt} = obj.integrators{iInt}.compute(nodalFunc);
            end
        end

        function f = integrateAndSum(obj,unfFun)
            f         = zeros(obj.dofs,1);
            iBoundary = 0;
            if (isequal(class(unfFun),'UnfittedBoundaryFunction'))
                bcMesh   = unfFun.unfittedMesh.boundaryCutMesh.mesh;
                obj.test = eval([class(obj.test),'.create(bcMesh,1)']);
            end
            for iInt = 1:obj.nInt
                integrator = obj.integrators{iInt};
                if contains(class(integrator),'Composite')
                    iBoundary = iBoundary + 1;
                    newUnfFun = unfFun.obtainFunctionAtExternalBoundary(iBoundary);
                    intLoc    = integrator.integrateAndSum(newUnfFun);
                    int       = obj.computeGlobalIntegralFromLocal(intLoc,iBoundary);
                elseif isequal(class(integrator), 'RHSintegrator_ShapeFunction')
                    int = integrator.compute(unfFun,obj.test);
                else
                    mesh    = obj.unfittedMesh.backgroundMesh;
                    testFun = eval([obj.testClass,'.create(mesh,1)']);
                    int     = integrator.compute(unfFun,testFun);
                end
                f = f + int;
            end
        end

    end

    methods (Access = private)

        function init(obj, cParams)
            obj.nInt         = numel(cParams.compositeParams);
            obj.unfittedMesh = cParams.unfittedMesh;
            obj.testClass    = class(cParams.test);
            mesh             = cParams.unfittedMesh.backgroundMesh;
            obj.test         = eval([obj.testClass,'.create(mesh,1)']);
            obj.dofs         = size(obj.test.fValues,1);
        end

        function createIntegrators(obj,cParams)
            params = cParams.compositeParams;
            for iInt = 1:obj.nInt
                s = params{iInt};
                s.quadType = cParams.quadType;
                s.test = cParams.test;
                integrator = RHSintegrator.create(s);
                obj.integrators{end+1} = integrator;
            end
        end

        function fun = createInnerFunction(obj, charFun)
            s.mesh    = obj.unfittedMesh.backgroundMesh;
            s.fValues = unique(charFun.evaluate(ones(2,1)))*ones(s.mesh.nnodes,1); % !!
            p1Old     = P1Function(s);
            innerMesh = obj.unfittedMesh.innerMesh.mesh;
            connecIG  = obj.unfittedMesh.innerMesh.globalConnec;
            fun       = p1Old.restrict2cell(innerMesh,connecIG);
        end

        function int = computeGlobalIntegralFromLocal(obj, intLoc, iBoundary)
            connecBG              = obj.unfittedMesh.unfittedBoundaryMesh.getGlobalConnec{iBoundary};
            meshes                = obj.unfittedMesh.unfittedBoundaryMesh.getActiveMesh();
            connecBL              = meshes{iBoundary}.backgroundMesh.connec;
            loc2glob(connecBL(:)) = connecBG(:);
            int                   = zeros(obj.dofs,1);
            int(loc2glob)         = intLoc;
        end

    end

end