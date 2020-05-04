classdef MeshAndGaussOfUnfittedQuadMesh < handle
    
    properties (Access = protected)
        backgroundMesh
        levelSet
    end
    
    methods (Access = protected)
        
        function init(obj)
            obj.computeBackgroundMesh();
            obj.computeLevelSet();
        end
        
        function computeLevelSet(obj)
            obj.levelSet = [-0.05 0.2 -0.5 0.1 0.1 -1 1 -0.2 -0.5]';
            %ls = [-0.05 0.2 -0.5 -0.1 0.1 -1 1 -0.2 -0.5]';
        end
        
        function computeBackgroundMesh(obj)
            s.coord = [0 0;1 0;1 1;0 1;2 0;2 1;0 2;1 2;2 2];
            s.connec = [1 2 3 4; 2 5 6 3; 4 3 8 7; 3 6 9 8];            
            m = Mesh().create(s);
            obj.backgroundMesh = m;
        end
        
        
        
        function [connecFull,connecCut,cutElems] = computeConnecCutAndFull(obj,ls,connec)
            lsInElem = obj.computeLevelSetInElem(ls,connec);
            isFull  = all(lsInElem<0,2);
            isEmpty = all(lsInElem>0,2);
            isCut = ~isFull & ~isEmpty;
            connecFull = connec(isFull,:);
            connecCut  = connec(isCut,:);
            cutElems = find(isCut);
        end
        
        
        
    end
    
    methods (Access = protected, Static)
        
        function m = computeMesh(connec,coord)
            s.connec = connec;
            s.coord = coord;
            m = Mesh().create(s);
        end
        
        function lsElem = computeLevelSetInElem(ls,connec)
            lsElem(:,1) = ls(connec(:,1));
            lsElem(:,2) = ls(connec(:,2));
            lsElem(:,3) = ls(connec(:,3));
            lsElem(:,4) = ls(connec(:,4));
        end
        
    end
    
    
    
    
    
    
    
end