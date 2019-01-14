classdef MeshPlotter_Interior_2D < MeshPlotter_Abstract
    properties (Access = private)
        meshUnfitted
        meshBackground
        fullCells
    end
    
    methods (Access = public)
        function plot(obj,meshUnfitted,ax)
            obj.init(meshUnfitted);
            obj.plotCutCells(ax);
            obj.plotInteriorCells(ax);
        end
        
        function meshUnfitted = patchRemovedDimension(obj,meshUnfitted,removedDim,removedDimCoord)
            meshUnfitted = obj.patchCutCellsCoords(meshUnfitted,removedDim,removedDimCoord);
            meshUnfitted = obj.patchFullCellsCoords(meshUnfitted,removedDim,removedDimCoord);
        end
    end
    
    methods (Access = private)
        function init(obj,meshUnfitted)
            obj.meshUnfitted = meshUnfitted;
            obj.meshBackground = meshUnfitted.meshBackground;
            obj.fullCells = meshUnfitted.background_full_cells;
        end
        
        function plotCutCells(obj,ax)
            obj.plotSurface(ax,obj.meshUnfitted.coord,obj.meshUnfitted.connec);
        end
        
        function plotInteriorCells(obj,ax)
            obj.plotSurface(ax,obj.meshBackground.coord,obj.meshBackground.connec(obj.fullCells,:));
        end
        
        function meshUnfitted = patchCutCellsCoords(obj,meshUnfitted,removedDim,removedDimCoord)
            newCoords = obj.patchCoords(meshUnfitted,removedDim,removedDimCoord);
            meshUnfitted.changeCoordinates(newCoords);
        end
        
        
        function meshUnfitted = patchFullCellsCoords(obj,meshUnfitted,removedDim,removedDimCoord)
            newCoords = obj.patchCoords(meshUnfitted.meshBackground,removedDim,removedDimCoord);
            meshUnfitted.meshBackground.changeCoordinates(newCoords);
        end
        
        function newCoords = patchCoords(obj,mesh,removedDim,removedDimCoord)
            nnode = size(mesh.coord,1);
            newCoords = zeros(nnode,3);
            newCoords = obj.addRemovedDimension(newCoords,removedDim,removedDimCoord,nnode);
            newCoords = obj.addPreservedDimensions(newCoords,mesh.coord,removedDim);
        end
    end
    
    methods (Access = private, Static)
        function plotSurface(ax,coord,connec)
            hold on
            patch(ax,'vertices',coord,'faces',connec,...
                'edgecolor',[0.5 0 0], 'edgealpha',0.5,'edgelighting','flat',...
                'facecolor',[1 0 0],'facelighting','flat')
        end
        
        function newCoords = addRemovedDimension(newCoords,removedDim,removedDimCoord,nnode)
            patch = removedDimCoord*ones(1,nnode);
            newCoords(:,removedDim) = patch;
        end
        
        function newCoords = addPreservedDimensions(newCoords,oldCoords,removedDim)
            iOld = 0;
            for iNew = 1:3
                if iNew ~= removedDim
                    iOld = iOld+1;
                    newCoords(:,iNew) = oldCoords(:,iOld);
                end
            end
        end
    end
end