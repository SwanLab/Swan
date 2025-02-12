classdef PatchedMeshPlotter_Abstract < handle
    
    properties (Access = protected)
        meshUnfitted
        meshBackground
        fullCells
    end
    
    methods (Access = protected, Static, Abstract)
        plotSurface(ax,coord,connec)
    end
    
    methods (Access = public)
        
        function plot(obj,meshUnfitted,ax,bF)
            obj.init(meshUnfitted,bF);
            obj.plotCutCells(ax);
            obj.plotInteriorCells(ax);
        end
        
        function meshUnfitted = patchRemovedDimension(obj,meshUnfitted,removedDim,removedDimCoord)
            obj.init(meshUnfitted);            
            meshUnfitted = obj.patchCutCellsCoords(meshUnfitted,removedDim,removedDimCoord);
            meshUnfitted = obj.patchFullCellsCoords(meshUnfitted,removedDim,removedDimCoord);
        end
    end
    
    methods (Access = private)
        
        function init(obj,meshUnfitted,backgroundFullCells)
            obj.meshUnfitted = meshUnfitted;
            obj.meshBackground = meshUnfitted.meshBackground;
            obj.fullCells = backgroundFullCells;
        end
        
        function plotCutCells(obj,ax)
            obj.plotSurface(ax,obj.meshUnfitted.coord,obj.meshUnfitted.connec);
        end
        
        function plotInteriorCells(obj,ax)
            obj.plotSurface(ax,obj.meshBackground.coord,obj.meshBackground.connec(obj.fullCells,:));
        end
        
        function meshUnfitted = patchCutCellsCoords(obj,meshUnfitted,removedDim,removedDimCoord)
            newCoords = obj.patchCoords(meshUnfitted,removedDim,removedDimCoord);
            meshUnfitted.setCoord(newCoords);
        end
        
        
        function meshUnfitted = patchFullCellsCoords(obj,meshUnfitted,removedDim,removedDimCoord)
            newCoords = obj.patchCoords(meshUnfitted.meshBackground,removedDim,removedDimCoord);
            meshUnfitted.meshBackground.setCoord(newCoords);
        end
        
        function newCoords = patchCoords(obj,mesh,removedDim,removedDimCoord)
            nnode = size(mesh.coord,1);
            newCoords = zeros(nnode,mesh.ndim+1);
            newCoords = obj.addRemovedDimension(newCoords,removedDim,removedDimCoord,nnode);
            newCoords = obj.addPreservedDimensions(newCoords,mesh.coord,removedDim);
        end
        
        function newCoords = addPreservedDimensions(obj,newCoords,oldCoords,removedDim)
            iOld = 0;
            for iNew = 1:obj.meshBackground.ndim+1
                if iNew ~= removedDim
                    iOld = iOld+1;
                    newCoords(:,iNew) = oldCoords(:,iOld);
                end
            end
        end
        
    end
    
    methods (Access = private, Static)
        
        function newCoords = addRemovedDimension(newCoords,removedDim,removedDimCoord,nnode)
            patch = removedDimCoord*ones(1,nnode);
            newCoords(:,removedDim) = patch;
        end
        
    end
end

