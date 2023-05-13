classdef LevelSetPeriodicAndOriented < LevelSetCreator

    properties (Access = private)
        epsilon
        epsilons
        dilation
        dilatedOrientation
        cellCoord
        phi
        y1
        y2
        deformedCoord
        m1
        m2
        meshD
    end

    properties (Access = private)
        mesh
        remesher
        cellLevelSetParams
        orientationVector
        nCells
    end

    methods (Access = public)

        function obj = LevelSetPeriodicAndOriented(cParams)
            obj.init(cParams);
            obj.computeDilation();
            obj.computeDilatedOrientationVector();
            obj.createDeformedCoord();
            obj.createRemesher();
            obj.interpolateDeformedCoord();
            obj.interpolateDilatation();
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
            obj.orientationVector  = cParams.orientationVector;
            obj.cellLevelSetParams = cParams.cellLevelSetParams;
        end

        function computeDilation(obj)
            s.orientationVector = obj.orientationVector;
            s.mesh  = obj.mesh;
            dC = DilationComputer(s);
            d  = dC.compute();
            obj.dilation = d;
        end

        function computeDilatedOrientationVector(obj)
            s.fValues = exp(obj.dilation.fValues);
            s.mesh    = obj.mesh;
            er = P1Function(s);
            for iDim = 1:obj.mesh.ndim
                b  = obj.orientationVector{iDim};
                dO = P1Function.times(er,b);
                obj.dilatedOrientation{iDim} = dO;
            end
        end        

        function createDeformedCoord(obj)
            s.mesh               = obj.mesh;
            s.orientationVector  = obj.orientationVector;
            s.dilatedOrientation = obj.dilatedOrientation;
            c = ConformalMappingComputer(s);
            defCoord = c.compute();
            obj.deformedCoord = defCoord.getFvaluesAsVector();
        end

        function createRemesher(obj)
            s.mesh    = obj.mesh.createDiscontinuousMesh();
            s.nLevels = 2;
            r  = Remesher(s);
            r.remesh();
            obj.remesher = r;
        end

        function createCellCoord(obj)
            nDim = obj.mesh.ndim;
            for iDim = 1:nDim
                x = obj.deformedCoord(:,iDim);
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
            nDim  = obj.mesh.ndim;
            nnode = obj.mesh.nnodeElem;
            for iDim = 1:nDim
                yD = obj.deformedCoord(:,iDim);
                y(1,:,:) = reshape(yD,nnode,[]);
                yI = obj.interpolateDiscontinousFunction(y);
                yI = abs(yI);
                yT(:,iDim) = yI;
            end
            obj.deformedCoord = yT;
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
            r  = obj.dilation.fValues;
            r  = obj.interpolateContinousFunctionToDisc(r);
            obj.dilation = r;
        end

        function t = computeMinLengthInUnitCell(obj)
            r = obj.dilation;
            hC = obj.epsilon*exp(-r);
            hmin = min(hC);
            hmax = max(hC);
            % hcut = (hmax+hmin)/0.6;%/4;%/2;
            hcut = 0;%0.00001*obj.epsilon;
            t = hcut./hC;
        end

        function fV = createDiscontinousValues(obj,r)
            s.mesh    = obj.mesh;
            s.fValues = r;
            fC = P1Function(s);
            fD = fC.project('P1D');
            fV = fD.fValues;
        end

        function vq = interpolateFunction(obj,v,mesh)
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

        function vq = interpolateDiscontinousFunction(obj,v)
            s.mesh    = obj.mesh;
            s.fValues = v;
            f         = P1DiscontinuousFunction(s);
            r = obj.remesher;
            f = r.interpolate(f);
            vq = f.getFvaluesAsVector();
        end

        function [v] = refine(obj,m,v)
            s.mesh    = m;
            s.fValues = v;
            f         = P1DiscontinuousFunction(s);
            r = obj.remesher;
            f = r.interpolate(f);
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
