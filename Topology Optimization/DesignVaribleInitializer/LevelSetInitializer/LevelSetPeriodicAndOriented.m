classdef LevelSetPeriodicAndOriented < LevelSetCreator
 
    properties (Access = private)
        epsilon
        epsilons
        dilation
        cellCoord
        phi
        y1
        y2
        m1
        m2
        meshD
    end

    properties (Access = private)
        mesh
        remesher
        cellLevelSetParams
        theta
        nCells
    end

    methods (Access = public)
        
        function obj = LevelSetPeriodicAndOriented(cParams)
            obj.init(cParams); 
            obj.meshD = obj.mesh.createDiscontinuousMesh();
            obj.computeDilation();            
            obj.createMapping();            
            obj.computeDeformedCoord();  
  
            obj.interpolateDilatation();            
            obj.interpolateM1M2();
        end

        function ls = computeLS(obj)
            nEps = length(obj.epsilons);
            ls = cell(nEps,1);
            for iEps = 1:nEps
                obj.epsilon = obj.epsilons(iEps);
                obj.computeLevelSet();
                ls{iEps} = obj.getValue();
             end                        
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
           obj.theta              = cParams.theta;           
           obj.epsilons             = cParams.epsilons;
           obj.cellLevelSetParams = cParams.cellLevelSetParams; 
        end   

        function computeDilation(obj)
            s.theta = obj.theta;
            s.mesh  = obj.mesh;
            dC = DilationFieldComputer(s);
            d  = dC.compute();
            %dC.plot();
            obj.dilation = d;
        end      

        function createMapping(obj)
            s.mesh     = obj.mesh;
            s.theta    = obj.theta;
            s.dilation = obj.dilation;
            c = ConformalMappingComputer(s);
            phiV = c.compute();
            % c.plot();
            obj.phi = phiV;
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
        
        function computeDeformedCoord(obj)
            y1(1,:,:) = reshape(obj.phi(:,1),3,[]);
            y2(1,:,:) = reshape(obj.phi(:,2),3,[]);
            y1 = obj.interpolateDiscontinousFunction(y1);
            y2 = obj.interpolateDiscontinousFunction(y2);
            obj.y1 = abs(y1);
            obj.y2 = abs(y2);            
        end 

        function interpolateM1M2(obj)
            m1 = obj.cellLevelSetParams.widthH;
            m2 = obj.cellLevelSetParams.widthV;
            obj.m1 = obj.interpolateContinousFunctionToDisc(m1);
            obj.m2 = obj.interpolateContinousFunctionToDisc(m2);
        end
        
        function thresholdParameters(obj)
            mL = obj.computeMinLengthInUnitCell();
            s.minLengthInUnitCell = mL;            
            t = MparameterThresholder(s);
            m1 = t.thresh(obj.m1);
            m2 = t.thresh(obj.m2);            
            obj.cellLevelSetParams.widthH = m1;
            obj.cellLevelSetParams.widthV = m2; 
        end

        function fDI = interpolateContinousFunctionToDisc(obj,fC)
            fD   = obj.createDiscontinousValues(fC);
            fDI  = obj.interpolateDiscontinousFunction(fD);   
        end

        function interpolateDilatation(obj)
            r  = obj.dilation;            
            r  = obj.interpolateContinousFunctionToDisc(r);   
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

        function fV = createDiscontinousValues(obj,r)
            m = obj.mesh;
            s.type    = m.type;
            s.connec  = m.connec;
            s.fValues = r;
            mD = obj.meshD;
            fC = P1Function(s);
            fD = fC.createP1Discontinous(mD);            
            fV = fD.fValues;
        end
        
        function vq = interpolateDiscontinousFunction(obj,v)
            m = obj.meshD;
            s.type    = m.type;
            s.connec  = m.connec;
            s.fValues = v;
            f         = P1DiscontinuousFunction(s);
            r = obj.remesher;
            f = r.interpolate(f);            
            vq = f.getFvaluesAsVector();
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
           f = abs(cos(pi/2*(y))).^2;
         %  f = y - floor(y);
        end                        
        
    end
    
end