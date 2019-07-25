classdef Mesh_Unfitted_Composite < Mesh_Unfitted
    
    properties (Access = public)
        meshBackground
    end
    
    properties (GetAccess = public, SetAccess = private)
        meshes
        nMeshes
        
        innerMeshOLD
        boxFaceMeshes
        globalConnectivities
        activeBoxFaceMeshesList
        activeMeshesList
        nActiveBoxFaces
        nActiveMeshes
        unfittedType = 'COMPOSITE'
        
        coord
        connec
        subcellIsoCoords
        cellContainingSubcell
    end
    
    properties (Access = private)
        nodesInBoxFaces
        isBoxFaceMeshActive
        nboxFaces
        ndim
        nsides = 2;
        
        totalMesh
        
        removedDimensions
        removedDimensionCoord
    end
    
    methods (Access = public)
        
        function obj = Mesh_Unfitted_Composite(cParams)
            obj.init(cParams.meshBackground);
            obj.createInteriorMesh(cParams);
            obj.createBoxMeshes();
            obj.createMeshes();
        end
        
        function computeMesh(obj,levelSet)
            obj.levelSet_background = levelSet;
            obj.computeInteriorMesh(levelSet);
            obj.computeBoxMeshes(levelSet);
            obj.createMeshes();
        end
        
        function M = computeMass(obj)
            Mi = obj.computeInteriorMass();
            Mb = obj.computeBoxMass();
            M = Mi + Mb;
        end
        
        function plot(obj)
            h = figure;
            obj.add2plot(axes(h));
            light
            axis equal off
            hold off
        end
        
        function add2plot(obj,ax)
            obj.innerMeshOLD.add2plot(ax);
            for iActive = 1:obj.nActiveBoxFaces
                iFace = obj.activeBoxFaceMeshesList(iActive);
                obj.boxFaceMeshes{iFace}.add2plot(ax,obj.removedDimensions(iFace),obj.removedDimensionCoord(iFace));
            end
        end
        
        function aMeshes = getActiveMeshes(obj)
            aMeshes = cell(1,obj.nActiveMeshes);
            for iActive = 1:obj.nActiveMeshes
                iMesh = obj.activeMeshesList(iActive);
                aMeshes{iActive} = obj.meshes{iMesh};
            end
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,meshBackground)
            obj.ndim = meshBackground.ndim;
            obj.meshBackground = meshBackground;
            obj.nboxFaces = obj.ndim*obj.nsides;
            
            s.coord = meshBackground.coord;
            s.connec = meshBackground.connec;
            obj.totalMesh = Mesh_Total(s);
            
            obj.removedDimensions = obj.totalMesh.removedDimensions;
            obj.removedDimensionCoord = obj.totalMesh.removedDimensionCoord;
        end
        
        function createMeshes(obj)
            obj.meshes{1} = obj.innerMeshOLD;
            for iface = 1:obj.nboxFaces
                obj.meshes{iface+1} = obj.boxFaceMeshes{iface};
            end
            obj.nMeshes = numel(obj.meshes);
        end
        
        function createInteriorMesh(obj,cParams)
            s.coord  = obj.totalMesh.coord;
            s.connec = obj.totalMesh.connec;
            s.meshBackground = obj.totalMesh.innerMeshOLD;
            s.interpolationBackground = Interpolation.create(s.meshBackground,'LINEAR');
            s.unfittedType = cParams.unfittedType;
            obj.innerMeshOLD = Mesh_Unfitted_Single(s);
        end
        
        function computeInteriorMesh(obj,levelSet)
            obj.innerMeshOLD.computeMesh(levelSet);
        end
        
        function M = computeInteriorMass(obj)
            M = obj.innerMeshOLD.computeMass();
        end
        
        function createBoxMeshes(obj)
            iFace = 0;
            fMeshes = obj.totalMesh.boxFaceMeshes;
            fNodes = obj.totalMesh.nodesInBoxFaces;
            fGlobalConnec = obj.totalMesh.globalConnectivities;
            for idime = 1:obj.ndim
                for iside = 1:obj.nsides
                    iFace = iFace + 1;
                    mesh = fMeshes{iFace};
                    nodesInBoxFace = fNodes{iFace};
                    obj.boxFaceMeshes{iFace}        = obj.createBoxFaceMesh(mesh);
                    obj.nodesInBoxFaces{iFace}      = nodesInBoxFace;
                    obj.globalConnectivities{iFace} = fGlobalConnec{iFace};
                end
            end
        end
        
        function computeBoxMeshes(obj,levelSet)
            obj.resetBoxMeshes();
            iface = 0;
            obj.isBoxFaceMeshActive = false([1 obj.nboxFaces]);
            for idime = 1:obj.ndim
                for iside = 1:obj.nsides
                    iface = iface + 1;
                    boxFaceMesh = obj.boxFaceMeshes{iface};
                    mshBack = boxFaceMesh.meshBackground;
                    lsBoxFace = levelSet(obj.nodesInBoxFaces{iface});
                    if obj.isBoxMeshActive(lsBoxFace,mshBack)
                        obj.boxFaceMeshes{iface}.compute(lsBoxFace);
                        obj.isBoxFaceMeshActive(iface) = true;
                    end
                end
            end
            obj.updateActiveBoxFaceMeshesList();
        end
        
        function resetBoxMeshes(obj)
            obj.createBoxMeshes();
        end
        
        function M = computeBoxMass(obj)
            M = 0;
            for iactive = 1:obj.nActiveBoxFaces
                iface = obj.activeBoxFaceMeshesList(iactive);
                M = M + obj.boxFaceMeshes{iface}.computeMass();
            end
        end
        
        function boxFaceMesh = createBoxFaceMesh(obj,mesh)
            interp = Interpolation.create(mesh,'LINEAR');
            s.unfittedType = 'INTERIOR';
            s.meshBackground = mesh;
            s.interpolationBackground = interp;
            cParams = SettingsMeshUnfitted(s);
            boxFaceMesh = UnfittedMesh(cParams);
        end
        
        
        function indexes = findConnecIndexes(obj,coord_indexes,nnode)
            number_of_valid_nodes_per_element = sum(ismember(obj.meshBackground.connec,coord_indexes),2);
            indexes = number_of_valid_nodes_per_element == nnode;
        end
        
        function updateActiveBoxFaceMeshesList(obj)
            obj.activeBoxFaceMeshesList = find(obj.isBoxFaceMeshActive);
            obj.activeMeshesList = find([true obj.isBoxFaceMeshActive]);
            obj.nActiveBoxFaces = length(obj.activeBoxFaceMeshesList);
            obj.nActiveMeshes = length(obj.activeMeshesList);
        end
    end
    
    methods (Static, Access = private)
        
        function itIs = isBoxMeshActive(levelSet,meshBack)
            phi_nodes = levelSet(meshBack.connec);
            phi_case = sum((sign(phi_nodes)<0),2);
            itIs = (any(phi_case));
        end
        
    end
    
    methods
        
        function coord = get.coord(obj)
            coord = obj.innerMeshOLD.coord;
        end
        
        function connec = get.connec(obj)
            connec = obj.innerMeshOLD.connec;
        end
        
        function s = get.subcellIsoCoords(obj)
            s = obj.innerMeshOLD.subcellIsoCoords;
         end
        
         function c = get.cellContainingSubcell(obj)
            c = obj.innerMeshOLD.cellContainingSubcell;
        end
        
    end
    
end

