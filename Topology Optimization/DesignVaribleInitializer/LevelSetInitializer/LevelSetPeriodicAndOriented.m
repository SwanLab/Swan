classdef LevelSetPeriodicAndOriented < handle

    properties (Access = private)
        epsilon
        epsilons
        cellCoord
        phi
        y1
        y2
        deformedCoord
        m1
        m2
    end

    properties (Access = private)
        mesh
        remesher
        cellLevelSetParams
        orientationVectors
        nCells
    end

    methods (Access = public)
        
        function obj = LevelSetPeriodicAndOriented(cParams)
            obj.init(cParams);
            obj.createDeformedCoord();
            obj.createRemesher();
            obj.interpolateDeformedCoord();
            obj.interpolateM1M2();
        end

        function ls = computeLS(obj,epsilons)
            nEps = length(epsilons);
            ls = cell(nEps,1);
            for iEps = 1:nEps
                obj.epsilon = epsilons(iEps);
                lsF = obj.computeLevelSet();
                ls{iEps} = lsF;
            end
        end

        function mF = getFineMesh(obj)
            fMesh = obj.remesher.fineMesh;
            mF = fMesh.createDiscontinuousMesh();
        end

    end

    methods (Access = protected)

        function ls = computeLevelSet(obj)
            obj.createCellCoord();
            obj.thresholdParameters();
            ls = obj.createCellLevelSet();
        end

    end


    methods (Access = private)

        function init(obj,cParams)
            obj.mesh               = cParams.mesh;
            obj.orientationVectors = cParams.orientationVectors;
            obj.cellLevelSetParams = cParams.cellLevelSetParams;
        end

        function createDeformedCoord(obj)
            c  = obj.orientationVectors;
            dC = c.computeDeformedCoordinates();
            obj.deformedCoord = dC;
        end

        function createRemesher(obj)
            s.mesh    = obj.mesh.createDiscontinuousMesh();
            s.nLevels =  3;
            r  = Remesher(s);
            r.remesh();
            obj.remesher = r;
        end

        function createCellCoord(obj)
            nDim = obj.mesh.ndim;
            for iDim = 1:nDim
                x = obj.deformedCoord.fValues(:,iDim);
                y = obj.computeMicroCoordinate(x);
                y = obj.periodicFunction(y);
                yT(:,iDim) = y;
            end
            obj.cellCoord = yT;
        end

        function ls = createCellLevelSet(obj)
            s       = obj.cellLevelSetParams;
            s.xCoorCenter = 0;%obj.cellCoord(:,1);
            s.yCoorCenter = 0;%obj.cellCoord(:,2);
            m = obj.getFineMesh();
            g = GeometricalFunction(s);
            ls = g.computeLevelSetFunction(m);
        end

        function interpolateDeformedCoord(obj)
            y = obj.deformedCoord;
            y = obj.interpolateDiscontinousFunction(y);
            y.fValues = abs(y.fValues);
            obj.deformedCoord = y;
        end

        function interpolateM1M2(obj)
            m1 = obj.cellLevelSetParams.xSide;
            m2 = obj.cellLevelSetParams.ySide;               
            obj.m1 = obj.interpolateContinousFunctionToDisc(m1);
            obj.m2 = obj.interpolateContinousFunctionToDisc(m2);
          %  p  = obj.cellLevelSetParams.pnorm;                 
          %  p = obj.interpolateContinousFunctionToDisc(p); 
          %  obj.cellLevelSetParams.pnorm = p;
        end

        function thresholdParameters(obj)
            mL = obj.computeMinLengthInUnitCell();
            s.minLengthInUnitCell = mL;
            t = MparameterThresholder(s);
            obj.m1.fValues = t.thresh(obj.m1.fValues);
            obj.m2.fValues = t.thresh(obj.m2.fValues);
            obj.cellLevelSetParams.xSide = @(x) obj.m1.evaluate(x);
            obj.cellLevelSetParams.ySide = @(x) obj.m2.evaluate(x);
        end

 
% 
%         function interpolateDilatation(obj)
%             r  = obj.orientationVectors.getDilatation();
%             r  = obj.interpolateContinousFunctionToDisc(r);
%             obj.dilation = r;
%         end

        function t = computeMinLengthInUnitCell(obj)
%            r = obj.dilation;
%            hC = obj.epsilon*exp(-r);
%            hmin = min(hC);
%            hmax = max(hC);
%             hcut = (hmax+hmin)/4;%/4;%/2;
%            %hcut = 0;%0.00001*obj.epsilon;
%            t = hcut./hC;
            t = 0;
        end

       function fDI = interpolateContinousFunctionToDisc(obj,fC)
            fD   = obj.createDiscontinousValues(fC);
            fDI  = obj.interpolateDiscontinousFunction(fD);
        end        

        function fV = createDiscontinousValues(obj,f)
            s.mesh    = obj.mesh;
            s.fValues = f;
            s.order   = 'P1';
            fC = LagrangianFunction(s);
            fD = fC.project('P1D');
            fV = fD;
        end

        function f = interpolateDiscontinousFunction(obj,v)
            f = v;
            r = obj.remesher;
            f = r.interpolate(f);
         %   vq = f.getFvaluesAsVector();
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
