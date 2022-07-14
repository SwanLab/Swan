classdef LevelSetPeriodicAndOriented < LevelSetCreator
 
    properties (Access = private)
        epsilon
        cellCoord
        mapping
    end
    
    properties (Access = private)
      mesh
      backgroundMesh
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
            obj.createMapping();
            obj.createCellCoord();
            obj.thresholdParameters();
            obj.createCellLevelSet();            
        end 
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
           obj.mesh               = cParams.mesh;
           obj.backgroundMesh     = cParams.backgroundMesh;
           obj.angle              = cParams.angle;
           obj.cellLevelSetParams = cParams.cellLevelSetParams; 
           obj.nCells             = cParams.nCells;
        end
        
        function createEpsilon(obj)
            L = obj.mesh.computeCharacteristicLength();
            obj.epsilon = L/obj.nCells;            
        end
        
        function createMapping(obj)
            s.mesh  = obj.mesh;
            s.theta = obj.angle;
            map = ConformalMappingComputer(s);
            map.compute();
            map.plot();            
            obj.mapping = map;
        end        
        
        function createCellCoord(obj)
            [y1,y2] = obj.applyMapping();                                                                                                                                                                  
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
        
        function [y1,y2] = applyMapping(obj)
            y1 = obj.mapping.phi(:,1);
            y2 = obj.mapping.phi(:,2);
            y1 = obj.interpolateFunction(y1);
            y2 = obj.interpolateFunction(y2);
        end 
        
        function thresholdParameters(obj)
            mL = obj.computeMinLengthInUnitCell();
            s.minLengthInUnitCell = mL;
            m1 = obj.cellLevelSetParams.widthH;
            m2 = obj.cellLevelSetParams.widthV;
            t = MparameterThresholder(s);
            m1 = t.thresh(m1);
            m2 = t.thresh(m2);            
            obj.cellLevelSetParams.widthH = m1;
            obj.cellLevelSetParams.widthV = m2;                        
        end
        
        function t = computeMinLengthInUnitCell(obj)
            r = obj.mapping.dilation;
            r = obj.interpolateFunction(r);            
            hC = obj.epsilon*exp(-r);
            hmin = min(hC);
            hmax = max(hC);
           % hcut = (hmax+hmin)/0.6;%/4;%/2;
            hcut = 0.00001*obj.epsilon;
            t = hcut./hC;                           
        end
        
        function vq = interpolateFunction(obj,v)
            X = obj.mesh.coord(:,1);
            Y = obj.mesh.coord(:,2);
            F = scatteredInterpolant(X,Y,v);
            xB = obj.backgroundMesh.coord(:,1);
            yB = obj.backgroundMesh.coord(:,2);
            vq = F(xB,yB);
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