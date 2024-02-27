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
            uMeshFun  = unfFun.unfittedMeshFunction;
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

                    % % %
                    ss.mesh     = obj.unfittedMesh.innerMesh.mesh;
                    ss.type     = 'ShapeFunction';
                    ss.quadType = 'LINEAR';
                    int2        = RHSintegrator.create(ss);
                    test2       = LagrangianFunction.create(ss.mesh, obj.test.ndimf, obj.test.order);
                    int2        = int2.compute(uMeshFun.innerMeshFunction,test2);
                    int2        = obj.computeGlobalIntegralFromLocalInner(int2);
                    % % %
                else
                    s.mesh      = obj.unfittedMesh.backgroundMesh;
                    s.feFunType = obj.testClass;
                    s.ndimf     = obj.test.ndimf;
                    testFun     = LagrangianFunction.create(s.mesh, s.ndimf, obj.test.order);
                    int         = integrator.compute(unfFun,testFun);

                    % % %
                    ss.mesh     = obj.unfittedMesh.innerCutMesh.mesh;
                    ss.type     = 'ShapeFunction';
                    ss.quadType = 'LINEAR';
                    int2        = RHSintegrator.create(ss);
                    test2       = LagrangianFunction.create(ss.mesh, obj.test.ndimf, obj.test.order);
                    int2        = int2.compute(uMeshFun.innerCutMeshFunction,test2);
                    %int2        = obj.computeGlobalIntegralFromLocalInnerCut(int2,uMeshFun.innerNodes,uMeshFun.cutEdges,uMeshFun.nonCutMesh);
                    % % %
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

        function int = computeGlobalIntegralFromLocalInner(obj, intLoc)
            connecBG              = obj.unfittedMesh.innerMesh.globalConnec;
            connecBL              = obj.unfittedMesh.innerMesh.mesh.connec;
            loc2glob(connecBL(:)) = connecBG(:);
            int                   = zeros(obj.dofs,1);
            int(loc2glob)         = intLoc;
        end

        function int = computeGlobalIntegralFromLocalInnerCut(obj, intLoc, innerNodes, cutEdges, cMeshGlobal)
            intLocInner           = intLoc(1:length(innerNodes));
            intLocInnerCut        = intLoc(length(innerNodes)+1:end);
            int                   = zeros(obj.dofs,1);
            int(innerNodes)       = intLocInner;
            interp = Interpolation.create('LINE','LINEAR');
            xIso = cutEdges.xIso';
            NxV  = interp.computeShapeFunctions(xIso)';
            nonCutNodes = unique(cMeshGlobal.connec(:));
            cutNodes1 = nonCutNodes(cutEdges.edges(:,1));
            cutNodes2 = nonCutNodes(cutEdges.edges(:,2));
            int(cutNodes1) = int(cutNodes1) + intLocInnerCut.*NxV(:,1);
            int(cutNodes2) = int(cutNodes2) + intLocInnerCut.*NxV(:,2);
        end

    end

end