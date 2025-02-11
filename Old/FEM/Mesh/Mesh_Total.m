classdef Mesh_Total < Mesh_Composite 
    
    properties (GetAccess = public, SetAccess = private)
        nBoxFaces
        
        innerMeshOLD
        boxFaceMeshes
        nodesInBoxFaces
        

        
        nnodes
        nnodeElem
        embeddedDim

    end
    
    properties (Access = private)
        borderNodes
        borderElements
        isExteriorMeshExplicit
    end
    
    methods (Access = public)
        
        function obj = Mesh_Total(cParams)
            obj.init(cParams);
            obj.createInteriorMesh();
            obj.createBoxFaceMeshes();
            obj.defineActiveMeshes();
            obj.type = obj.innerMeshOLD.type;
            obj.nelem = size(obj.connec,1);
            obj.nnodes = obj.innerMeshOLD.nnodes;
            obj.nnodeElem = obj.innerMeshOLD.nnodeElem;
            obj.createInterpolation();
            obj.computeElementCoordinates();
        end
        
        function S = computeMeanCellSize(obj)
            S = obj.innerMeshOLD.computeMeanCellSize();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.coord  = cParams.coord;
            obj.connec = cParams.connec;
            obj.obtainExteriorMesh(cParams);
            obj.ndim   = size(obj.coord,2);
            obj.embeddedDim = obj.ndim;
        end
        
        function obtainExteriorMesh(obj,cParams)
            obj.isExteriorMeshExplicit = false;
            if isfield(cParams,'borderNodes')
                if ~isempty(cParams.borderNodes)
                    obj.borderNodes    = cParams.borderNodes;
                    obj.borderElements = cParams.borderElements;
                    obj.isExteriorMeshExplicit = true;
                end
            end
        end
        
        function createInteriorMesh(obj)
            s.connec = obj.connec;
            s.coord  = obj.coord;
            obj.innerMeshOLD = Mesh.create(s);
            obj.append(obj.innerMeshOLD);
        end
        
        function createBoxFaceMeshes(obj)
            if obj.isExteriorMeshExplicit
                obj.computeExteriorMeshesFromData();
            else
               obj.computeExteriorMeshesFromBoxSides();
            end
        end
        
        function computeExteriorMeshesFromData(obj)
            s.backgroundMesh = obj.innerMeshOLD;
            s.borderNodes    = obj.borderNodes;
            s.borderElements = obj.borderElements;
            bC = BoundaryMeshCreatorFromRectangularBox(s);
            bMeshes = bC.create();
            obj.nBoxFaces = numel(bMeshes);
            for iM = 1:obj.nBoxFaces
                m = bMeshes{iM};
                obj.boxFaceMeshes{iM} = m;
                obj.append(m);
            end
        end
        
        function computeExteriorMeshesFromBoxSides(obj)
            s.backgroundMesh = obj.innerMeshOLD;
            s.dimension = 1:s.backgroundMesh.ndim;
            s.type = 'FromReactangularBox';
            bC = BoundaryMeshCreator.create(s);
            bMeshes = bC.create();
            obj.nBoxFaces = numel(bMeshes);
            for iM = 1:obj.nBoxFaces
                m = bMeshes{iM};
                obj.boxFaceMeshes{iM} = m;
                obj.append(m);
            end
        end
        
        function defineActiveMeshes(obj)
            obj.activeMeshesList = find([false true(1,obj.nBoxFaces)]);
            obj.nActiveMeshes     = numel(obj.activeMeshesList);
        end
        
    end
    
end