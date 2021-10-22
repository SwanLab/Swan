classdef UnfittedMeshSplitter < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        
    end
    
    properties (Access = private)
        uMesh
    end
    
    methods (Access = public)
        
        function obj = UnfittedMeshSplitter(cParams)
            obj.init(cParams)
        end
        
        function split(obj)
            innerMeshes    = obj.computeInnerMeshComponents();
            innerCutMeshes = obj.computeInnerCutMeshComponents();
            
            
            allMeshes = [innerMeshes, innerCutMeshes];
            remainingMeshes = allMeshes;
            
            iMeshA = 1;
            iComp = 1;
            while length(remainingMeshes)> 1
                compMesh = remainingMeshes(1);
                remainingMeshes = remainingMeshes(2:end);
                t = 0;
                while t< length(compMesh)
                meshA = compMesh(t+1);
                meshA = meshA{1};    
                switch class(meshA)
                    case 'InnerMesh'
                        nodesA = meshA.globalConnec(:);
                    case 'InnerCutMesh'
                        nodesA = meshA.mesh.connec(:);
                end
                k = 1;
                sameComp = [];
                for iMeshB = 1:length(remainingMeshes)
                    meshB = remainingMeshes(iMeshB);
                    meshB = meshB{1};    
                    
                    switch class(meshB)
                        case 'InnerMesh'
                            nodesB = meshB.globalConnec(:);
                        case 'InnerCutMesh'
                            nodesB = meshB.mesh.connec(:);
                    end
                    if obj.areSharingNodes(nodesA,nodesB)
                        sameComp(k) = iMeshB;
                        k = k+1;
                        compMesh{end+1}= meshB;
                    end
                end
               % sameComp = [1,sameComp];
                isRem = setdiff(1:length(remainingMeshes),sameComp);
                
               % compMesh = remainingMeshes(sameComp);
                remainingMeshes = remainingMeshes(isRem);
                t = t+1;
                end
                a{iComp} = compMesh;
                iComp = iComp+1;
                
            end
            
            %for
            
            
            
            
            
            for i = 1:length(a)
                aP = a(i);
                aP = aP{1};
                figure()
                hold on
                for j = 1:length(aP)
                    mi = aP{j};
                    s.edgeAlpha = 0;
                    s.faceAlpha = 1;
                    s.faceColor = 'black';
                    s.mesh = mi.mesh;
                    s = SettingsMeshPlotter(s);
                    
                    mP = MeshPlotter(s);
                    mP.plot();
                   % axis([0 2 -1 1])
                    %axis([-1 1 0 1])
                   % mi.mesh.plot();
                end
            end
            
            
            
            
            
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.uMesh = cParams.unfittedMesh;
        end
        
        function mI = computeInnerMeshComponents(obj)
            inner = obj.uMesh.innerMesh;
            faces = inner.mesh.connec;
            s.faces = faces;
            sp = SplitterInConnectedComponents(s);
            CC = sp.split();
            isComponent = CC(faces(:,1))';
            nComp = unique(isComponent);
            for iComp = 1:size(nComp)
                isComp = isComponent == nComp(iComp);
                
                s.fullCells = inner.fullCells(isComp);
                s.backgroundMesh = inner.backgroundMesh;
                mI{iComp} = InnerMesh(s);
                
                
            end
        end
        
        function mI = computeInnerCutMeshComponents(obj)
            innerCut   = obj.uMesh.innerCutMesh;
            faces = innerCut.mesh.connec;
            s.faces = faces;
            sp = SplitterInConnectedComponents(s);
            CC = sp.split();
            isComponent = CC(faces(:,1))';
            nComp = unique(isComponent);
            for iComp = 1:size(nComp)
                isComp = isComponent == nComp(iComp);
                
                s.connec = faces(isComp,:);
                s.coord   = innerCut.mesh.coord;
                s.mesh                  = Mesh(s);
                s.xCoordsIso            = innerCut.xCoordsIso(:,:,isComp);
                s.cellContainingSubcell = innerCut.cellContainingSubcell(isComp,:);
                mI{iComp} = InnerCutMesh(s);
                
                
                
            end
        end
        
    end
    
    methods (Static, Access = private)
        
        function theyAre = areSharingNodes(nodesA,nodesB)
            nodesAB = intersect(nodesA,nodesB);
            theyAre = ~isempty(nodesAB);
        end
        
    end
    
end