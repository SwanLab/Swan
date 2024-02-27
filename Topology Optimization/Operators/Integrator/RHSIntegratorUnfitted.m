classdef RHSIntegratorUnfitted < handle

    properties (Access = private)
        quadType
        unfittedMesh
    end

    properties (Access = private)
        innerIntegrator
   end

    methods (Access = public)
        function obj = RHSIntegratorUnfitted(cParams)
            obj.init(cParams);
            obj.createIntegratorInner();
        end

        function int = compute(obj,uMeshFun,test)
            fInner   = uMeshFun.innerMeshFunction;
            intInner = obj.integrateInnerMesh(fInner,test);
            int  = intInner+0;
        end
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.quadType     = cParams.quadType;
            obj.unfittedMesh = cParams.mesh;
        end

        function createIntegratorInner(obj)
            s.mesh     = obj.unfittedMesh.innerMesh.mesh;
            s.type     = 'ShapeFunction';
            s.quadType = obj.quadType;
            int        = RHSintegrator.create(s);
            obj.innerIntegrator = int;
        end

        function int = integrateInnerMesh(obj,f,test)
            intLoc        = obj.integrateInnerMeshLocal(f,test);
            iMesh         = obj.unfittedMesh.innerMesh;
            dofs          = size(test.fValues,1);
            int           = zeros(dofs,1);
            conG          = iMesh.globalConnec;
            conBL         = iMesh.mesh.connec;
            l2g(conBL(:)) = conG(:);
            int(l2g)      = intLoc;
        end

        function intLoc = integrateInnerMeshLocal(obj,f,test)
            m       = obj.unfittedMesh.innerMesh.mesh;
            testLoc = LagrangianFunction.create(m,test.ndimf,test.order);
            intLoc  = obj.innerIntegrator.compute(f,testLoc);
        end

    end
end