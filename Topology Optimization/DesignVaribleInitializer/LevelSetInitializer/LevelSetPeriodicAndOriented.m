classdef LevelSetPeriodicAndOriented < handle

    
    properties (Access = private)
        epsilon
        epsilons
        phi
        y1
        y2
        deformedCoord
        m1
        m2
        fineMesh
    end

    properties (Access = private)
        mesh
        remesher
        orientationVectors
        nCells
        nRemeshLevels
    end

    methods (Access = public)
        
        function obj = LevelSetPeriodicAndOriented(cParams)
            obj.init(cParams);
            obj.createDeformedCoord();
            obj.remeshAndInterpolate();
        end

        function ls = computeLS(obj,epsilons)
            nEps = length(epsilons);
            ls = cell(nEps,1);
            for iEps = 1:nEps
                eps = epsilons(iEps);
                lsF = obj.computeLevelSet(eps);
                ls{iEps} = lsF;
            end
        end

        function fM = getFineMesh(obj)
           fM = obj.fineMesh;
        end

    end

    methods (Access = protected)

        function ls = computeLevelSet(obj,eps)
            obj.thresholdParameters();
            ls = obj.createCellLevelSet(eps);
        end

    end


    methods (Access = private)

        function init(obj,cParams)
            obj.mesh               = cParams.mesh;
            obj.orientationVectors = cParams.orientationVectors;
            obj.nRemeshLevels      = cParams.nRemeshLevels;
            obj.m1                 = cParams.m1;
            obj.m2                 = cParams.m2;
        end

        function createDeformedCoord(obj)
            c  = obj.orientationVectors;
            dC = c.computeDeformedCoordinates();
            obj.deformedCoord = dC;
        end

        function remeshAndInterpolate(obj)
            m = obj.mesh;
            for iLevel = 1:obj.nRemeshLevels
                m = m.remesh();
                obj.deformedCoord = obj.deformedCoord.refine(m);
                obj.m1            = obj.m1.refine(m);
                obj.m2            = obj.m2.refine(m);                
            end
            obj.fineMesh = m;
        end

        function yT = evaluateCellCoord(obj,xV,eps)
            nDim = obj.mesh.ndim;
            xD = obj.deformedCoord.evaluate(xV);            
            for iDim = 1:nDim
                x = squeezeParticular(xD(iDim,:,:),1);
                y = obj.computeMicroCoordinate(x,eps);
                y = obj.periodicFunction(y);
                yT(iDim,:,:) = y;
            end
        end

        function ls = createCellLevelSet(obj,eps)
            s.mesh    = obj.fineMesh;
            s.fHandle = @(xV) obj.rectangle(xV,eps);
            s.ndimf    = 1;
            f  = AnalyticalFunction(s);
            ls = f.project('P1');
        end

        function fH = rectangle(obj,xV,eps)
            sx = obj.m1.evaluate(xV);
            sy = obj.m2.evaluate(xV);
            x0 = eps;
            y0 = eps;
            x = obj.evaluateCellCoord(xV,eps);
            x1 = x(1,:,:);
            x2 = x(2,:,:);
            p = 32;
         %   fH = ((abs(x1-x0)./sx).^p+(abs(x2-y0)./sy).^p).^(1/p) - 0.5;
            fH = max(abs(x1-x0)./sx,abs(x2-y0)./sy) - 0.5;
            fH = -fH;
        end

        function thresholdParameters(obj)
            mL = obj.computeMinLengthInUnitCell();
            s.minLengthInUnitCell = mL;
            t = MparameterThresholder(s);
            obj.m1.fValues = t.thresh(obj.m1.fValues);
            obj.m2.fValues = t.thresh(obj.m2.fValues);
        end

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

        function [y1,y2] = transformToFastCoord(obj,x1,x2)
            y1 = obj.computeMicroCoordinate(x1);
            y2 = obj.computeMicroCoordinate(x2);
        end

        function y = computeMicroCoordinate(obj,x,eps)
            y = (x-min(x(:)))/eps ;
        end

        function  [y1,y2] = makeCoordPeriodic(obj,y1,y2)
            y1 = obj.periodicFunction(y1);
            y2 = obj.periodicFunction(y2);
        end

    end

    methods (Access = private, Static)

        function f = periodicFunction(y)
            f = abs(cos(pi*(y)));%.^2;           
         %   f = (y - floor(y));
        end

    end

end
