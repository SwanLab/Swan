classdef BoundayConditionsPlotter_Abstract < handle
    
    properties (GetAccess = protected, SetAccess = private)
        axes
        mesh
        bc
        scaleD
        scaleN
        
        iD
        iN
        vD
        vN
        
        colorD = [0 .5 1]
        colorN = [.2 .8 .2]
        lineWidth = 2
        maxHeadSize = .5
    end
    
    methods (Access = public, Abstract)
        
        plotDirichlet(obj)
        plotNeumann(obj)
        
    end
    
    methods (Access = public)
        
        function plot(obj)
            obj.plotDirichlet();
            obj.plotNeumann();
        end
        
    end
    
    methods (Access = protected)
        
        function init(obj,axes,mesh,bc)
            obj.axes = axes;
            obj.mesh = mesh;
            obj.bc   = bc;
            
            obj.getBoundaryConditions();
            
            obj.computeScaleFactor();
        end
        
    end
    
    methods (Access = private)
        
        function getBoundaryConditions(obj)
            obj.findNodes();
            obj.findValues();
        end
        
        function findNodes(obj)
            obj.iN = obj.getNodesFromInputData(obj.bc.pointload);
            obj.iD = obj.getNodesFromInputData(obj.bc.dirichlet);
        end
        
        function findValues(obj)
            obj.vN = obj.getValuesFromInputData(obj.bc.pointload);
            obj.vD = obj.getValuesFromInputData(obj.bc.dirichlet);
        end
        
        function computeScaleFactor(obj)
            f = 0.2;
            
            maxCoord = max(obj.mesh.coord(:));
            minCoord = min(obj.mesh.coord(:));
            
            maxDirValue = max(abs(obj.vD(:)));
            maxNeuValue = max(abs(obj.vN(:)));
            
            obj.scaleD = f*(maxCoord-minCoord)/maxDirValue;
            obj.scaleN = f*(maxCoord-minCoord)/maxNeuValue;
        end
        
    end
    
    
    methods (Access = private, Static)
        
        function i = getNodesFromInputData(d)
            if ~isempty(d)
                i = unique(d(:,1));
            else
                i = [];
            end
        end
        
        function v = getValuesFromInputData(d)
            v = zeros(size(d));
            iArray = [];
            for i = 1:size(d,1)
                if ~isempty(iArray)
                    isNew = ~any(d(i,1) == iArray);
                else
                    isNew = true;
                end
                if isNew
                    iArray(end+1) = d(i,1);
                end
                index = find(iArray,d(i,1));
                v(index,d(i,2)) = d(i,3);
            end
            v(length(iArray)+1:end,:) = [];
        end
        
    end
    
end