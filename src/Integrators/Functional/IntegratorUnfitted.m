classdef IntegratorUnfitted < handle

    properties (Access = private)
        quadType
        unfittedMesh
    end

    properties (Access = private)
        innerIntegrator
        innerCutIntegrator
        boundCutIntegrator
        unfittedBoundIntegrator
    end

    methods (Access = public)
        function obj = IntegratorUnfitted(cParams)
            obj.init(cParams);
            obj.createIntegratorInner();
            obj.createIntegratorInnerCut();
            obj.createIntegratorBoundaryCut();
            obj.createIntegratorUnfittedBoundary();
        end

        function int = compute(obj,uFun)
            switch class(uFun)
                case 'UnfittedFunction'
                    intInner    = obj.integrateInnerMeshFunction(uFun);
                    intInnerCut = obj.integrateInnerCutMeshFunction(uFun);
                    int         = intInner+intInnerCut;
                case 'UnfittedBoundaryFunction'
                    intBoundCut = obj.integrateBoundaryCutMeshFunction(uFun);
                    intUnfBound = obj.integrateUnfittedBoundaryMeshFunction(uFun);
                    int         = intBoundCut+intUnfBound;
            end
        end
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.quadType     = cParams.quadType;
            obj.unfittedMesh = cParams.mesh;
        end

        function createIntegratorInner(obj)
            if ~isempty(obj.unfittedMesh.innerMesh)
                m   = obj.unfittedMesh.innerMesh.mesh;
                int = Integrator.create('Function',m,obj.quadType);
                obj.innerIntegrator = int;
            end
        end

        function createIntegratorInnerCut(obj)
            if ~isempty(obj.unfittedMesh.innerCutMesh)
                m   = obj.unfittedMesh.innerCutMesh.mesh;
                int = Integrator.create('Function',m,obj.quadType);
                obj.innerCutIntegrator = int;
            end
        end

        function createIntegratorBoundaryCut(obj)
            if ~isempty(obj.unfittedMesh.boundaryCutMesh)
                m   = obj.unfittedMesh.boundaryCutMesh.mesh;
                int = Integrator.create('Function',m,obj.quadType);
                obj.boundCutIntegrator = int;
            end
        end

        function createIntegratorUnfittedBoundary(obj)
            int = @(uMesh) Integrator.create('Unfitted',uMesh,obj.quadType);
            obj.unfittedBoundIntegrator = int;
        end

        function int = integrateInnerMeshFunction(obj,uFun)
            if ~isempty(obj.unfittedMesh.innerMesh)
                fInner = uFun.innerMeshFunction;
                int    = obj.innerIntegrator.compute(fInner);
            else
                int = 0;
            end
        end

        function int = integrateInnerCutMeshFunction(obj,uFun)
            if ~isempty(obj.unfittedMesh.innerCutMesh)
                fInnerCut = uFun.innerCutMeshFunction;
                int       = obj.innerCutIntegrator.compute(fInnerCut);
            else
                int = 0;
            end
        end

        function int = integrateBoundaryCutMeshFunction(obj,uFun)
            if ~isempty(obj.unfittedMesh.boundaryCutMesh)
                fBoundCut = uFun.boundaryCutMeshFunction;
                int       = obj.boundCutIntegrator.compute(fBoundCut);
            else
                int = 0;
            end
        end

        function int = integrateUnfittedBoundaryMeshFunction(obj,uFun)
            if ~isempty(obj.unfittedMesh.unfittedBoundaryMesh)
                uFunBound  = uFun.unfittedBoundaryMeshFunction.activeFuns;
                uMeshBound = obj.unfittedMesh.unfittedBoundaryMesh.getActiveMesh();
                nFun = length(uFunBound);
                int  = 0;
                for i = 1:nFun
                    funi       = uFunBound{i}.backgroundFunction;
                    uMeshi     = uMeshBound{i};
                    s.fun      = funi;
                    s.uMesh    = uMeshi;
                    uF         = UnfittedFunction(s);
                    integrator = obj.unfittedBoundIntegrator(uMeshi);
                    inti       = integrator.compute(uF);
                    int        = int + inti;
                end
            else
                int = 0;
            end
        end
    end
end