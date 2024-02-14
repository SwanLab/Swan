classdef LevelSetPeriodicAndOriented < LevelSetCreator

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
%            obj.interpolateDilatation();
            obj.interpolateM1M2();
        end

        function ls = computeLS(obj,epsilons)
            nEps = length(epsilons);
            ls = cell(nEps,1);
            for iEps = 1:nEps
                obj.epsilon = epsilons(iEps);
                obj.computeLevelSet();
                ls{iEps} = obj.getValue();
            end
        end

        function mF = getFineMesh(obj)
            fMesh = obj.remesher.fineMesh;
            mF = fMesh.createDiscontinuousMesh();
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

        function createCellLevelSet(obj)
            s       = obj.cellLevelSetParams;
            s.coord = obj.cellCoord;
            ls = LevelSetCreator.create(s);
            obj.levelSet = ls.getValue();
        end

        function interpolateDeformedCoord(obj)
            y = obj.deformedCoord;
            y = obj.interpolateDiscontinousFunction(y);
            y = abs(y);
            obj.deformedCoord.fValues = y;
        end

        function interpolateM1M2(obj)
            m1 = obj.cellLevelSetParams.widthH;
            m2 = obj.cellLevelSetParams.widthV;               
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
            m1 = t.thresh(obj.m1);
            m2 = t.thresh(obj.m2);
            obj.cellLevelSetParams.widthH = m1;
            obj.cellLevelSetParams.widthV = m2;
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

        function vq = interpolateDiscontinousFunction(obj,v)
            f = v;
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
