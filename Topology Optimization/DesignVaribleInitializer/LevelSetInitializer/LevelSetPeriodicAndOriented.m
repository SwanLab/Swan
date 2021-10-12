classdef LevelSetPeriodicAndOriented < LevelSetCreator
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        epsilon
        cellCoord
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
            obj.createCellCoord();
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
            s.mesh  = obj.mesh;
            s.theta = obj.angle;
            map = ConformalMappingComputer(s);
            map.compute();
            map.plot();
            y1 = map.phi(:,1);
            y2 = map.phi(:,2);
            
            
          %  r = map.dilation;
            r = obj.interpolateFunction(map.dilation);
            
            hC = obj.epsilon*exp(-r);
            
            m1 = obj.cellLevelSetParams.widthH;
            m2 = obj.cellLevelSetParams.widthV;
            
            hmin = min(hC);
            hmax = max(hC);
            hcut = (hmax+hmin)/4;%/2;
     %       hcut = 0.05*obj.epsilon;
  
            
           % m1 = obj.thresh(hcut,hC,m1);
           % m2 = obj.thresh(hcut,hC,m2);
            
            obj.cellLevelSetParams.widthH = m1;
            obj.cellLevelSetParams.widthV = m2;
            
            
            y1 = obj.interpolateFunction(y1);
            y2 = obj.interpolateFunction(y2);
            
        end   
        
        function m = thresh(obj,hcut,hC,m)     
            
            alpha = hcut./hC;
    
            isTooSmall = alpha < 2; 
            
            p = sum(isTooSmall)/size(isTooSmall,1);
            
            isMSmall05 = m < 0.5;
            isMLarge05 = ~isMSmall05;
            
            
            m(isTooSmall & isMSmall05) = 0;
            m(isTooSmall & isMLarge05) = 1;                    

            isHcLargerThan2Hmin = alpha < 1/2;
            isMSmallerThan1_HminHc = m < alpha;
            isMLargerThan2_HminHc = m > 0.5*alpha;
            itIs = isHcLargerThan2Hmin & isMSmallerThan1_HminHc & isMLargerThan2_HminHc;
            m(itIs) = alpha(itIs);
            
            isMLargerThan2_HminHc = m < 0.5*alpha;
            itIs = isHcLargerThan2Hmin & isMLargerThan2_HminHc;
            m(itIs) = 0;         
            
            
            isMSmallerThan1_HminHc = m > 1- alpha;
            isMLargerThan2_HminHc = m < 1 - 0.5*alpha;
            itIs = isHcLargerThan2Hmin & isMSmallerThan1_HminHc & isMLargerThan2_HminHc;
            m(itIs) = 1- alpha(itIs);
            
            isMLargerThan2_HminHc = m > 1 - 0.5*alpha;
            itIs = isHcLargerThan2Hmin & isMLargerThan2_HminHc;
            m(itIs) = 1;         
            
            
            
            
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