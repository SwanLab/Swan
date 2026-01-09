classdef UnfittedMeshSplitter < handle
    
    properties (Access = public)
        meshGroups
    end
    
    properties (Access = private)
        allMeshes
        meshGroupLabels
    end
    
    properties (Access = private)
        uMesh
    end
    
    methods (Access = public)
        
        function obj = UnfittedMeshSplitter(cParams)
            obj.init(cParams)
        end
        
        function split(obj)
            obj.computeAllMeshes();
            obj.computeMeshGroupLabels();
            obj.computeMeshesGroups();
        end
        
        function plot(obj)
            m = obj.meshGroups;
            for i = 1:length(m)
                mG = m{i};
                obj.plotMeshGroup(mG);
            end
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.uMesh = cParams.unfittedMesh;
        end
        
        function computeAllMeshes(obj)
            innerMeshes    = obj.computeInnerMeshes();
            innerCutMeshes = obj.computeInnerCutMeshes();
            aM = [innerMeshes; innerCutMeshes];
            obj.allMeshes = aM;
        end
        
        function m = computeInnerMeshes(obj)
            mI = obj.uMesh.innerMesh;
            m = obj.computeSubMeshComponents(mI);                      
        end
        
        function m = computeInnerCutMeshes(obj)
            mI = obj.uMesh.innerCutMesh;
            m = obj.computeSubMeshComponents(mI);                      
        end        
        
        function computeMeshesGroups(obj)
            mLabel = obj.meshGroupLabels;
            nG = max(mLabel);
            aM = obj.allMeshes;
            m = cell(nG,1);
            for iG = 1:nG
                isG = iG == mLabel;
                m{iG} = aM(isG);
            end
            obj.meshGroups = m;
        end
        
        function computeMeshGroupLabels(obj)
            aM = obj.allMeshes;
            A = obj.computeMeshSharingNodesMatrix(aM);
            d = graph(A,'omitselfloops');
            gL = conncomp(d);
            obj.meshGroupLabels = gL;
        end
        
        function A = computeMeshSharingNodesMatrix(obj,meshes)
            nM = length(meshes);
            A = false(nM,nM);
            for iM = 1:nM
                mI = meshes(iM);
                nodI = obj.obtainNodes(mI);
                for jM = (iM+1):nM
                    mJ = meshes(jM);
                    nodJ = obj.obtainNodes(mJ);
                    isSharing = obj.areSharingNodes(nodI,nodJ);
                    A(iM,jM) = isSharing;
                end
            end
            A = A+A';
        end

        function mI = computeSubMeshComponents(obj,subMesh)
            m         = subMesh;
            faces     = m.mesh.connec;
            compLabel = obj.computeComponentLabel(faces);
            nComp = unique(compLabel);
            mI    = cell(length(nComp),1);
            for iComp = 1:length(nComp)
                isComp = compLabel == nComp(iComp);
                sM = obj.createSubMeshComponent(m,isComp);
                mI{iComp} = sM;
            end
        end    

        function plotMeshGroup(obj,mG)
            figure()
            hold on
            for j = 1:length(mG)
                s.mesh = mG{j}.mesh;
                s.edgeAlpha = 0;
                s.faceAlpha = 1;
                s.faceColor = 'black';
                s.isBackground = false;
                mP = MeshPlotter(s);
                mP.plot();
                xmax = max(obj.uMesh.backgroundMesh.coord(:,1));
                xmin = min(obj.uMesh.backgroundMesh.coord(:,1));
                ymax = max(obj.uMesh.backgroundMesh.coord(:,2));
                ymin = min(obj.uMesh.backgroundMesh.coord(:,2));
                axis([xmin xmax ymin ymax])
            end

        end
        

    end
    
    methods (Static, Access = private)
        
        function sM = createSubMeshComponent(m,nodes)
            switch class(m)
                case 'InnerMesh'
                    s.fullCells      = m.fullCells(nodes);
                    s.backgroundMesh = m.backgroundMesh;
                    sM = InnerMesh(s);
                case 'InnerCutMesh'
                    s.connec  = m.mesh.connec(nodes,:);
                    s.coord   = m.mesh.coord;
                    s.mesh                  = Mesh.create(s);
                    s.xCoordsIso            = m.xCoordsIso(:,:,nodes);
                    s.cellContainingSubcell = m.cellContainingSubcell(nodes,:);
                    sM = InnerCutMesh(s);
            end

        end

       function compID = computeComponentLabel(faces)
            s.faces = faces;
            sp = SplitterInConnectedComponents(s);
            mComp = sp.split();
            compID = mComp(faces(:,1))';
        end            
        
        function nodes = obtainNodes(mesh)
            m = mesh{1};
            switch class(m)
                case 'InnerMesh'
                    nodes = m.globalConnec(:);
                case 'InnerCutMesh'
                    nodes = m.mesh.connec(:);
            end
        end
        
        function theyAre = areSharingNodes(nodesA,nodesB)
            nodesAB = intersect(nodesA,nodesB);
            theyAre = ~isempty(nodesAB);
        end
        
        
    end
    
end