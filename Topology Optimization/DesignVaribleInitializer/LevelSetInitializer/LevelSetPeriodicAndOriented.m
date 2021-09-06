classdef LevelSetPeriodicAndOriented < LevelSetCreator
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        epsilon
        cellCoord
    end
    
    properties (Access = private)
      coord
      angle        
      cellLevelSetParams
      nCells
    end
    
    methods (Access = public)
        
        function obj = LevelSetPeriodicAndOriented(cParams)
            obj.init(cParams); 
            obj.computeLevelSet();
        end
        
    end
    
    methods (Access = protected)
        
        function computeLevelSet(obj)
            obj.createEpsilon();            
            obj.createCellCoord();
            obj.createCellLevelSet();            
        end 
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
           obj.coord              = cParams.coord;
           obj.angle              = cParams.angle;
           obj.cellLevelSetParams = cParams.cellLevelSetParams; 
           obj.nCells             = cParams.nCells;
        end

        function createEpsilon(obj)
            xy = obj.coord;
            x1max = max(xy(:,1));
            x1min = min(xy(:,1));
            x2max = max(xy(:,2));
            x2min = min(xy(:,2));            
            xMax = max(x1max,x2max);
            xMin = min(x1min,x2min);
            obj.epsilon = (xMax-xMin)/obj.nCells;            
        end        
        
        function createCellCoord(obj)
            xy = obj.coord;            
            x1 = xy(:,1);
            x2 = xy(:,2);        
            y1 = x1;
            y2 = x2;
            [y1,y2] = obj.rotateCoordinates(y1,y2);                                                                        
            [y1,y2] = obj.transformToFastCoord(y1,y2);                                       
            [y1,y2] = obj.makeCoordPeriodic(y1,y2);                                   
            obj.cellCoord = [y1,y2];
        end        
        
        function createCellLevelSet(obj)
           s       = obj.cellLevelSetParams;
           s.coord = obj.cellCoord;
           ls = LevelSetCreator.create(s);
           obj.levelSet = ls.getValue;
        end
        
        function [y1,y2] = rotateCoordinates(obj,x1,x2)
            theta = obj.angle; 
            R(1,1,:) = cos(theta);
            R(1,2,:) = -sin(theta);
            R(2,1,:) = sin(theta);
            R(2,2,:) = cos(theta);
            y1 = (squeeze(R(1,1,:)).*x1 + squeeze(R(1,2,:)).*x2);
            y2 = (squeeze(R(2,1,:)).*x1 + squeeze(R(2,2,:)).*x2);           
        end    
        
        function [y1,y2] = transformToFastCoord(obj,x1,x2)
            y1 = obj.computeMicroCoordinate(x1);
            y2 = obj.computeMicroCoordinate(x2);
        end   
        
        function y = computeMicroCoordinate(obj,x)
            eps = obj.epsilon;        
            y = (x-min(x))/eps;
        end          
        
        function  [y1,y2] = makeCoordPeriodic(obj,y1,y2)
            y1 = obj.periodicFunction(y1);
            y2 = obj.periodicFunction(y2);
        end        
        
    end
    
    methods (Access = private, Static)
        
        function f = periodicFunction(y)
            %f = abs(cos(pi/2*y)).^2;
            f = y - floor(y);
        end                        
        
    end
    
end