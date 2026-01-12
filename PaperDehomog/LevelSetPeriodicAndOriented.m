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
        density
        r
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
            obj.r  = project(obj.orientationVectors.getDilation(),'P1D'); 
            obj.remeshAndInterpolate();
        end

        function ls = computeLS(obj,epsilons)
            nEps = size(epsilons,1);
            ls = cell(nEps,1);
            for iEps = 1:nEps
                eps = epsilons(iEps,:);
                lsF = obj.computeLevelSet(eps);
                ls{iEps} = lsF;
            end
        end

        function fM = getFineMesh(obj)
           fM = obj.fineMesh;
        end

        function rho = getFineDensity(obj)
           rho = obj.density;
        end

    end

    methods (Access = protected)

        function ls = computeLevelSet(obj,eps)
            obj.thresholdParameters(eps);
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
            obj.density            = cParams.density;
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
                obj.r             = obj.r.refine(m);  
                obj.density       = obj.density.refine(m);  
            end
            obj.fineMesh = m;
        end

        function y = evaluateCellCoord(obj,xV,eps)
            x = obj.deformedCoord.getFvaluesByElem();
            
          %  x = obj.deformedCoord.evaluate(xV);             
            y = obj.computeMicroCoordinate(x,eps);
            y = obj.periodicFunction(y);
        end

        function ls = createCellLevelSet(obj,eps)
     %       s.operation  = @(xV) obj.geometricalFunction(xV,eps);
    %        s.ndimf      = 1;
    %        s.mesh       = obj.fineMesh;
    %        f  = DomainFunction(s);
    %        ls = project(f,'P1');      

            x = obj.mesh.coordElem;

            fValues = obj.geometricalFunction(x,eps);
            s.mesh    = obj.fineMesh;
            s.fValues = fValues(:);
            s.order   = 'P1D';              
            ls = LagrangianFunction(s);  
            ls = project(ls,'P1');
        end

        function fH = geometricalFunction(obj,xV,eps)
            sx = obj.m1.getFvaluesByElem();
            sy = obj.m2.getFvaluesByElem();
            %sx = obj.m1.evaluate(xV);
            %sy = obj.m2.evaluate(xV);
            x  = obj.evaluateCellCoord(xV,eps);            
            s.xSide = sx;
            s.ySide = sy;
            s.xCoorCenter = 0;
            s.yCoorCenter = 0;
            s.pnorm = 4;
            s.type = 'SmoothRectangleInclusion';
            %s.type = 'RectangleInclusion';
            g = GeometricalFunction(s);
            f = g.getHandle;
            fH = f(x);
        end

        function thresholdParameters(obj,eps)
            mL = obj.computeMinLengthInUnitCell(eps);
            s.minLengthInUnitCell = mL;
            t = MparameterThresholder(s);
            obj.m1.setFValues(t.thresh(obj.m1.fValues));
            obj.m2.setFValues(t.thresh(obj.m2.fValues));
        end

        function t = computeMinLengthInUnitCell(obj,epsilon)
             hC = epsilon*exp(-obj.r);
%            hmin = min(hC);
%            hmax = max(hC);
%             hcut = (hmax+hmin)/4;%/4;%/2;
            hcut = 0.05*epsilon;
            t = hcut./hC;
            t = project(t,'P1D');
            t = t.fValues;
        end

        function y = computeMicroCoordinate(obj,x,eps)
            nDim = size(x,1);
            y = zeros(size(x));
            for iDim = 1:nDim
                xI    = x(iDim,:,:);
                %xImin = min(xI(:))
               % xImin = eps(iDim)
               % xImin = 0;
                y(iDim,:,:) = (xI)/(eps);
            end
            %y = (x-min(x(:))-eps)/eps;
        end

    end

    methods (Access = private, Static)

        function f = periodicFunction(y)
            %f = ((cos(1*pi*(y)))).^2;           
        %    f = (cos(2*pi*y));   

            f = abs(y-floor(y)-0.5);
        end

    end

end
