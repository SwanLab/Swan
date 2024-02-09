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
                bcMesh      = unfFun.unfittedMesh.boundaryCutMesh.mesh;
                s.mesh      = bcMesh;
                s.feFunType = obj.testClass;
                s.ndimf     = obj.test.ndimf;
                obj.test    = LagrangianFunction.create(s.mesh, s.ndimf, obj.test.order);
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
                    s.mesh      = obj.unfittedMesh.backgroundMesh;
                    s.feFunType = obj.testClass;
                    s.ndimf     = obj.test.ndimf;
                    testFun     = LagrangianFunction.create(s.mesh, s.ndimf, obj.test.order);
                    int         = integrator.compute(unfFun,testFun);
                end
                f = f + int;
            end
        end

    end

    methods (Access = private)

        function init(obj, cParams)
            obj.nInt         = numel(cParams.compositeParams);
            obj.unfittedMesh = cParams.unfittedMesh;
            obj.testClass    = cParams.test.order;
            s.mesh           = cParams.unfittedMesh.backgroundMesh;
            s.feFunType      = obj.testClass;
            s.ndimf          = cParams.test.ndimf;
            obj.test         = LagrangianFunction.create(s.mesh, s.ndimf, cParams.test.order);
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