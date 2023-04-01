classdef RHSintegrator_Composite < handle

    properties (GetAccess = public, SetAccess = private)
        integrators
        nInt
        npnod
    end

    properties (Access = private)
        RHScells
        RHSsubcells
        unfittedMesh
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

        function f = integrateAndSum(obj,nodalFunc)
            f = 0;
            for iInt = 1:obj.nInt
                integrator = obj.integrators{iInt};
                if contains(class(integrator),'Composite')
                    int = integrator.integrateAndSum(nodalFunc);
                elseif isequal(class(integrator), 'RHSintegrator_ShapeFunctionFun')
                    p1 = obj.createInnerP1(nodalFunc);
                    intLoc = integrator.compute(p1);
                    int = obj.computeGlobalIntegralFromLocal(intLoc);
                else
                    int = integrator.compute(nodalFunc);
                end
                f = f + int;
            end
        end

    end

    methods (Access = private)

        function init(obj, cParams)
            obj.nInt = numel(cParams.compositeParams);
            obj.npnod = cParams.npnod;
            obj.unfittedMesh = cParams.unfittedMesh;
        end

        function createIntegrators(obj,cParams)
            params = cParams.compositeParams;
            for iInt = 1:obj.nInt
                s = params{iInt};
                integrator = RHSintegrator.create(s);
                obj.integrators{end+1} = integrator;
            end
        end

        function p1 = createInnerP1(obj, F)
            innerMesh = obj.unfittedMesh.innerMesh;
            connecIG = innerMesh.globalConnec;
            connecIL  = innerMesh.mesh.connec;
            innerL2G(connecIL(:)) = connecIG(:);

            innerDofsGlobal = unique(connecIG(:));
            innerDofs = unique(connecIL);

            fV_global = zeros(length(F),1);
            fV_global(innerDofsGlobal) = F(innerDofsGlobal);

            fV_localInner = zeros(length(innerDofs),1);
            fV_localInner(innerDofs) = fV_global(innerL2G(innerDofs));

            s.mesh = innerMesh.mesh;
            s.fValues = fV_localInner;
            p1  = P1Function(s);
        end

        function int = computeGlobalIntegralFromLocal(obj, intLoc)
            innerMesh = obj.unfittedMesh.innerMesh;
            connecIG = innerMesh.globalConnec;
            connecIL  = innerMesh.mesh.connec;
            innerL2G(connecIL(:)) = connecIG(:);
            innerDofs = unique(connecIL);

            int = zeros(obj.npnod,1);
            int(innerL2G(innerDofs)) = intLoc(innerDofs);
        end

    end

end