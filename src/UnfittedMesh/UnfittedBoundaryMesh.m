classdef UnfittedBoundaryMesh < handle
    
    properties (Access = public)
        meshes
    end
    
    properties (Access = private)
        nBoundaries
        levelSet
        activeMeshes
        globalConnec
        nodesInBoxFaces
    end
    
    properties (Access = private)
        boundaryMesh
    end
    
    methods (Access = public)
        
        function obj = UnfittedBoundaryMesh(cParams)
            obj.init(cParams)
            obj.createUnfittedMeshes();
            obj.obtainGlobalConnec();
            obj.obtainNodesInBoxFaces();
        end
        
        function compute(obj,ls)
            obj.levelSet = ls;
            obj.computeActiveMesh();
            obj.computeUnfittedMeshes();
        end

        function uF = computeBoundaryMeshFunction(obj,fP1)
            uF = [];
            if ~isempty(obj.meshes)
                aMeshes = obj.getActiveMesh();
                for i = 1:length(aMeshes)
                    uMeshi     = aMeshes{i};
                    connecLoc  = uMeshi.backgroundMesh.connec;
                    if ~isempty(connecLoc)
                        connecGlob = obj.getGlobalConnec{i};
                        glob2loc(connecLoc(:)) = connecGlob(:);
                        s.fValues  = fP1.fValues(glob2loc);
                        s.mesh     = uMeshi.backgroundMesh;
                        s.order    = 'P1';
                        fbackMeshi = LagrangianFunction(s);
                       % glob2loc   = [];
                        fi         = uMeshi.obtainFunctionAtUnfittedMesh(fbackMeshi);
                        uF.activeFuns{i} = fi;
                    end
                end
            end
        end

        function m = getActiveMesh(obj)
            m = obj.getActiveField(obj.meshes);
        end
        
        function g = getGlobalConnec(obj)
            g = obj.getActiveField(obj.globalConnec);
        end
        
        function v = computeVolume(obj)
            m = obj.getActiveMesh();
            v = 0;
            for imesh = 1:numel(m)
                v = v + m{imesh}.computeMass();
            end
        end
        
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.boundaryMesh = cParams.boundaryMesh;
            obj.nBoundaries  = numel(obj.boundaryMesh);
            obj.activeMeshes = false(obj.nBoundaries,1);
        end
        
        function createUnfittedMeshes(obj)
            for iBoundary = 1:obj.nBoundaries
                bMesh = obj.boundaryMesh{iBoundary};
                s.backgroundMesh = bMesh.mesh;
                if ~isequal(bMesh.mesh.geometryType,'Line')
                    s.boundaryMesh = obj.createBoundaryMesh(bMesh);
                else
                    s.boundaryMesh = [];
                end
                obj.meshes{iBoundary} = UnfittedMesh(s);
            end
        end
        
        function obtainGlobalConnec(obj)
            for iBoundary = 1:obj.nBoundaries
                m = obj.boundaryMesh{iBoundary};
                obj.globalConnec{iBoundary} = m.globalConnec;
            end
        end
        
        function obtainNodesInBoxFaces(obj)
            for iBoundary = 1:obj.nBoundaries
                m = obj.boundaryMesh{iBoundary};
                obj.nodesInBoxFaces{iBoundary} = m.nodesInBoxFaces;
            end
        end
        
        function computeActiveMesh(obj)
            for iBoundary = 1:obj.nBoundaries
                isActive = obj.isUnfittedMeshActive(iBoundary);
                obj.activeMeshes(iBoundary) = isActive;
            end
        end
        
        function computeUnfittedMeshes(obj)
            for iBoundary = 1:obj.nBoundaries       
                isMeshActive = obj.activeMeshes(iBoundary);
                if isMeshActive
                   obj.computeUnfittedMesh(iBoundary);
                end
            end
        end
        
        function computeUnfittedMesh(obj,iBoundary)
            nodes = obj.nodesInBoxFaces{iBoundary};
            ls    = obj.levelSet(nodes);
            obj.meshes{iBoundary}.compute(ls);
        end
        
        function itIs = isUnfittedMeshActive(obj,iBoundary)
            nodes = obj.nodesInBoxFaces{iBoundary};
            ls    = obj.levelSet(nodes);
            itIs = any(sign(ls)<0);
        end
        
        function activeFields = getActiveField(obj,fields)
            activeFields = cell(0);
            iMesh = 1;
            for iBoundary = 1:obj.nBoundaries
                isMeshActive = obj.activeMeshes(iBoundary);
                if isMeshActive
                   activeFields{iMesh} = fields{iBoundary};
                   iMesh = iMesh +1;
                end
            end
        end
        
    end
    
    methods (Access = private, Static)
    
        function m = createBoundaryMesh(bMesh)
            s.backgroundMesh = bMesh.mesh;
            s.dimension      = setdiff(1:bMesh.mesh.ndim,bMesh.dimension);
            bC = BoundaryMeshCreatorFromRectangularBox(s);
            m = bC.create();
        end
        
    end
    
end