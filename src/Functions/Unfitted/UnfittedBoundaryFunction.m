classdef UnfittedBoundaryFunction < BaseFunction

    properties (Access = public)
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
            res = copy(obj1);
            switch class(v)
                case 'Test'
                    res.fun = DP(obj1.fun,v);
                    res.unfittedBoundaryMeshFunction = obj1.computeAtUnfittedBoundaryMesh(v);
                    f       = obj1.boundaryCutMeshFunction;
                    isoMesh = obj1.obtainIsoparametricMesh();
                    xV      = @(xVLoc) isoMesh.evaluate(xVLoc);
                    Ni      = DomainFunction.create(@(xVLoc) v.evaluate(xV(xVLoc)),f.mesh,1);
                    res.boundaryCutMeshFunction = DP(f,Ni);
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

        function uBMFun = computeAtUnfittedBoundaryMesh(obj,v)
            u      = obj.unfittedBoundaryMeshFunction.activeFuns;
            uBMFun = cell(size(u));
            for i = 1:length(u)
                uMeshBound = obj.unfittedMesh.unfittedBoundaryMesh.getActiveMesh();
                s.fun      = u{i}.backgroundFunction;
                s.uMesh    = uMeshBound{i};
                uF         = UnfittedFunction(s);
                uBMFun{i}  = DP(uF,v.updateMesh(s.uMesh.backgroundMesh));
            end
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

    methods (Access = protected)

        function evaluateNew(obj,xV)

        end

    end
end