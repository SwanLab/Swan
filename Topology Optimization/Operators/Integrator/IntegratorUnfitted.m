classdef IntegratorUnfitted < handle

    properties (Access = private)
        quadType
        unfittedMesh
    end

    properties (Access = private)
        innerIntegrator
        innerCutIntegrator
    end

    methods (Access = public)
        function obj = IntegratorUnfitted(cParams)
            obj.init(cParams);
            obj.createIntegratorInner();
            obj.createIntegratorInnerCut();
        end

        function int = compute(obj,uMeshFun)
            fInner      = uMeshFun.innerMeshFunction;
            fInnerCut   = uMeshFun.innerCutMeshFunction;
            intInner    = obj.innerIntegrator.compute(fInner);
            intInnerCut = obj.innerCutIntegrator.compute(fInnerCut);
            int         = intInner+intInnerCut;
        end
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.quadType     = cParams.quadType;
            obj.unfittedMesh = cParams.mesh;
        end

        function createIntegratorInner(obj)
            m   = obj.unfittedMesh.innerMesh.mesh;
            int = Integrator.create('Function',m,obj.quadType);
            obj.innerIntegrator = int;
        end

        function createIntegratorInnerCut(obj)
            m   = obj.unfittedMesh.innerCutMesh.mesh;
            int = Integrator.create('Function',m,obj.quadType);
            obj.innerCutIntegrator = int;
        end
    end
end