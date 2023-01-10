classdef LevelSetPeriodicAndOriented < LevelSetCreator
 
    properties (Access = private)
        epsilon
        dilation
        cellCoord
        phi
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
            obj.computeDilation();
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
        
        function computeDilation(obj)
            s.theta = obj.angle;
            s.mesh  = obj.mesh;
            dC = DilationFieldComputer(s);
            d  = dC.compute();
            %dC.plot();
            obj.dilation = d;
        end           
        
        function createMapping(obj)
            s.mesh     = obj.mesh;
            s.theta    = obj.angle;
            s.dilation = obj.dilation;
            c = ConformalMappingComputer(s);
            phiV = c.compute();
           % c.plot();            
            obj.phi = phiV;
        end        
        
        function createCellCoord(obj)
            [y1,y2] = obj.applyMapping();  
            y1 = abs(y1);
            y2 = abs(y2);
            [y1,y2] = obj.transformToFastCoord(y1,y2);                                           
            [y1,y2] = obj.makeCoordPeriodic(y1,y2);                                                          
            obj.cellCoord = [y1,y2];
        end        
        
        function createCellLevelSet(obj)
           s       = obj.cellLevelSetParams;
           s.coord = obj.cellCoord;
           ls = LevelSetCreator.create(s);
           obj.levelSet = ls.getValue();
        end
        
        function [y1,y2] = applyMapping(obj)
            y1(1,:,:) = reshape(obj.phi(:,1),3,[]);
            y2(1,:,:) = reshape(obj.phi(:,2),3,[]);
            mD = obj.mesh.createDiscontinuousMesh;
            y1 = obj.interpolateFunction(y1,mD);
            y2 = obj.interpolateFunction(y2,mD);
        end 
        
        function thresholdParameters(obj)
            mL = obj.computeMinLengthInUnitCell();
            s.minLengthInUnitCell = mL;
            m1 = obj.cellLevelSetParams.widthH;
            m2 = obj.cellLevelSetParams.widthV;
            mD = obj.mesh.createDiscontinuousMesh();
            m1 = obj.createDiscontinousValues(m1,obj.mesh,mD);
            m2 = obj.createDiscontinousValues(m2,obj.mesh,mD);
            m1 = obj.interpolateFunction(m1,obj.mesh);
            m2 = obj.interpolateFunction(m2,obj.mesh);
            t = MparameterThresholder(s);
            m1 = t.thresh(m1);
            m2 = t.thresh(m2);            
            obj.cellLevelSetParams.widthH = m1;
            obj.cellLevelSetParams.widthV = m2; 
       %     p = obj.cellLevelSetParams.pnorm;
       %     p = obj.createDiscontinousValues(p,obj.mesh,mD);            
       %     p = obj.interpolateFunction(p,obj.mesh);
       %     obj.cellLevelSetParams.pnorm = p;
        end
        
        function t = computeMinLengthInUnitCell(obj)
            r = obj.dilation;            
            mD = obj.mesh.createDiscontinuousMesh();
            rD = obj.createDiscontinousValues(r,obj.mesh,mD);
            r  = obj.interpolateFunction(rD,mD);            
            hC = obj.epsilon*exp(-r);
            hmin = min(hC);
            hmax = max(hC);
           % hcut = (hmax+hmin)/0.6;%/4;%/2;
            hcut = 0;%0.000001*obj.epsilon;
            t = hcut./hC;                           
        end

        function fV = createDiscontinousValues(obj,r,m,mD)
            s.type    = m.type;
            s.connec  = m.connec;
            s.fValues = r;
            fC = P1Function(s);
            fD = fC.createP1Discontinous(mD);            
            fV = fD.fValues;
        end
        
        function vq = interpolateFunction(obj,v,mesh)
%             m = mesh;
%             X = m.coord(:,1);
%             Y = m.coord(:,2);
%             F = scatteredInterpolant(X,Y,v);
%             xB = obj.backgroundMesh.coord(:,1);
%             yB = obj.backgroundMesh.coord(:,2);
%             vq = F(xB,yB);
           vq = obj.refine(mesh,v);
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

        function [v] = refine(obj,m,v)
            s.type    = m.type;
            s.connec  = m.connec;
            s.fValues = v;
            f         = P1DiscontinuousFunction(s);
            for i = 1:3
                mF = m.remesh();
                f  = f.refine(m,mF);
                m  = mF.createDiscontinuousMesh();
            end
            v = f.getFvaluesAsVector();
        end
        
    end
    
    methods (Access = private, Static)
        
        function f = periodicFunction(y)
           f = abs(cos(pi/2*(y))).^2;
         %  f = y - floor(y);
        end                        
        
    end
    
end