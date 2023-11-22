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
            f         = zeros(obj.dofs);
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
                    int       = integrator.integrateAndSum(newUnfFun);
                elseif isequal(class(integrator), 'RHSintegrator_ShapeFunction')
                    int = integrator.compute(unfFun,obj.test);
                else
                    int = integrator.compute(unfFun);
                end
                f = f + int;
            end
        end

    end

    methods (Access = private)

        function init(obj, cParams)
            obj.nInt         = numel(cParams.compositeParams);
            %obj.dofs         = cParams.npnod;
            obj.dofs         = size(cParams.test.fValues);
            obj.unfittedMesh = cParams.unfittedMesh;
            obj.testClass    = class(cParams.test); % EL HACK
            obj.test         = cParams.test;
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

        function int = computeGlobalIntegralFromLocal(obj, intLoc)
            innerMesh = obj.unfittedMesh.innerMesh;
            connecIG = innerMesh.globalConnec;
            connecIL  = innerMesh.mesh.connec;
            innerL2G(connecIL(:)) = connecIG(:);
            innerDofs = unique(connecIL);

            int = zeros(obj.dofs,1);
            int(innerL2G(innerDofs)) = intLoc(innerDofs);
        end

    end

end