classdef LevelSetPeriodicAndOriented < LevelSetCreator
 
    properties (Access = private)
        epsilon
        dilation
        cellCoord
        phi
        y1
        y2
        m1
        m2
    end

    properties (Access = private)
        mesh
        remesher
        backgroundMesh
        angle
        cellLevelSetParams
    end

    methods (Access = public)
        
        function obj = LevelSetPeriodicAndOriented(cParams)
            obj.init(cParams); 
            [y1,y2] = obj.applyMapping();  
            obj.y1 = abs(y1);
            obj.y2 = abs(y2);  
            obj.interpolateDilatation();
            obj.interpolateM1M2();
        end

        function computeLs(obj,epsilon)
           obj.epsilon = epsilon;            
           obj.computeLevelSet();
        end
        
    end
    
    methods (Access = protected)
        
        function computeLevelSet(obj)
            obj.createCellCoord();
            obj.thresholdParameters();
            obj.createCellLevelSet();
        end

    end
        
    
    methods (Access = private)
        
        function init(obj,cParams)
           obj.mesh               = cParams.mesh;
           obj.remesher           = cParams.remesher;           
           obj.backgroundMesh     = cParams.backgroundMesh;
           obj.dilation           = cParams.dilation;           
           obj.phi                = cParams.phi;
           obj.cellLevelSetParams = cParams.cellLevelSetParams; 
        end
                                                
        function createCellCoord(obj)            
            [y1,y2] = obj.transformToFastCoord(obj.y1,obj.y2);                                           
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

        function interpolateM1M2(obj)
            m1 = obj.cellLevelSetParams.widthH;
            m2 = obj.cellLevelSetParams.widthV;
            mD = obj.mesh.createDiscontinuousMesh();
            m1 = obj.createDiscontinousValues(m1,obj.mesh,mD);
            m2 = obj.createDiscontinousValues(m2,obj.mesh,mD);
            obj.m1 = obj.interpolateFunction(m1,obj.mesh);
            obj.m2 = obj.interpolateFunction(m2,obj.mesh);
        end
        
        function thresholdParameters(obj)
            mL = obj.computeMinLengthInUnitCell();
            s.minLengthInUnitCell = mL;            
            t = MparameterThresholder(s);
            m1 = t.thresh(obj.m1);
            m2 = t.thresh(obj.m2);            
            obj.cellLevelSetParams.widthH = m1;
            obj.cellLevelSetParams.widthV = m2; 
       %     p = obj.cellLevelSetParams.pnorm;
       %     p = obj.createDiscontinousValues(p,obj.mesh,mD);            
       %     p = obj.interpolateFunction(p,obj.mesh);
       %     obj.cellLevelSetParams.pnorm = p;
        end

        function interpolateDilatation(obj)
            r = obj.dilation;            
            mD = obj.mesh.createDiscontinuousMesh();
            rD = obj.createDiscontinousValues(r,obj.mesh,mD);
            r  = obj.interpolateFunction(rD,mD);   
            obj.dilation = r;
        end
        
        function t = computeMinLengthInUnitCell(obj)
            r = obj.dilation;
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
%             for i = 1:2
%                mF = m.remesh();
%                f  = f.refine(m,mF);
%                m  = mF.createDiscontinuousMesh();
%             end
            r = obj.remesher;
            f = r.interpolate(f);

%             mF = obj.fineMesh;
%             f  = f.refine(m0,mF);
            
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