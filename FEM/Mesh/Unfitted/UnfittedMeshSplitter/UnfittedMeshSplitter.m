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

    end
    
    methods (Static, Access = private)
        
        function sM = computeSubMeshComponents(subMesh)
            s.subMesh = subMesh;
            uS = SubUnfittedMeshSplitter.create(s);
            sM = uS.split();
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
        
        function plotMeshGroup(mG)
            figure()
            hold on
            for j = 1:length(mG)
                s.mesh = mG{j}.mesh;
                s.edgeAlpha = 0;
                s.faceAlpha = 1;
                s.faceColor = 'black';
                s = SettingsMeshPlotter(s);
                mP = MeshPlotter(s);
                mP.plot();
            end
        end
        
    end
    
end