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

        function y = evaluateCellCoord(obj,xV,eps)
            x = obj.deformedCoord.evaluate(xV);             
            y = obj.computeMicroCoordinate(x,eps);
            y = obj.periodicFunction(y);
        end

        function ls = createCellLevelSet(obj,eps)
            s.operation  = @(xV) obj.geometricalFunction(xV,eps);
            s.ndimf      = 1;
            f  = DomainFunction(s);
            ls = f.project('P1',obj.fineMesh);            
        end

        function fH = geometricalFunction(obj,xV,eps)
            sx = obj.m1.evaluate(xV);
            sy = obj.m2.evaluate(xV);
            x0 = 0.5;
            y0 = 0.5;
            x  = obj.evaluateCellCoord(xV,eps);            
            s.xSide = sx;
            s.ySide = sy;
            s.xCoorCenter = x0;
            s.yCoorCenter = y0;
            s.pnorm = 32;
            s.type = 'SmoothRectangleInclusion';
            g = GeometricalFunction(s);
            f = g.getHandle;
            fH = f(x);
        end

        function fH = rectangle(obj,xV,eps)
            sx = obj.m1.evaluate(xV);
            sy = obj.m2.evaluate(xV);
            x0 = 0.5;
            y0 = 0.5;
            x  = obj.evaluateCellCoord(xV,eps);
            x1 = x(1,:,:);
            x2 = x(2,:,:);
            p = 2;
            fH = ((abs(x1-x0)./(0.5*sx)).^p+(abs(x2-y0)./(0.5*sy)).^p).^(1/p) - 1;
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

        function y = computeMicroCoordinate(obj,x,eps)
            nDim = size(x,1);
            y = zeros(size(x));
            for iDim = 1:nDim
                xI    = x(iDim,:,:);
                xImin = min(xI(:));
                y(iDim,:,:) = (xI-xImin)/(eps);
            end
            %y = (x-min(x(:))-eps)/eps;
        end

    end

    methods (Access = private, Static)

        function f = periodicFunction(y)
         %   f = abs(cos(1*pi*(y))).^2;           
            f = (y - floor(y));
        end

    end

end
