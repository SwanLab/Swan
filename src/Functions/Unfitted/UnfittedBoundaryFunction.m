classdef UnfittedBoundaryFunction < handle

    properties (Access = public)
        ndimf
        unfittedMesh
        boundaryCutMeshFunction
        unfittedBoundaryMeshFunction
    end

    properties (Access = private)
        fun
    end

    methods (Access = public)
        function obj = UnfittedBoundaryFunction(cParams)
            obj.init(cParams);
            obj.computeUnfittedMeshFunction();
        end

        function res = DP(obj1,v)
            switch class(v)
                case 'Test'
                    f       = obj1.boundaryCutMeshFunction;
                    isoMesh = obj1.obtainIsoparametricMesh();
                    xV      = @(xVLoc) isoMesh.evaluate(xVLoc);
                    Ni      = DomainFunction.create(@(xVLoc) v.evaluate(xV(xVLoc)),f.mesh,1);
                    res     = DP(f,Ni);
            end
        end
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.unfittedMesh = cParams.uMesh;
            obj.fun          = cParams.fun;
            obj.ndimf        = cParams.fun.ndimf;
        end

        function computeUnfittedMeshFunction(obj)
            uMeshFun = obj.unfittedMesh.obtainFunctionAtUnfittedMesh(obj.fun);
            obj.boundaryCutMeshFunction = uMeshFun.boundaryCutMeshFunction;
            obj.unfittedBoundaryMeshFunction = uMeshFun.unfittedBoundaryMeshFunction;
        end

        function m = obtainIsoparametricMesh(obj)
            coord      = obj.unfittedMesh.boundaryCutMesh.xCoordsIso;
            nDim       = size(coord,1);
            nNode      = size(coord,2);
            nElem      = size(coord,3);
            msh.connec = reshape(1:nElem*nNode,nNode,nElem)';
            msh.type   = obj.unfittedMesh.boundaryCutMesh.mesh.type;
            s.fValues  = reshape(coord,nDim,[])';
            s.mesh     = msh;
            s.order    = 'P1';
            m          = LagrangianFunction(s);
        end
    end
end