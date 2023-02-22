classdef ConformalMappingComputer < handle

    properties (Access = private)
        phi
        interpolator
        phiMapping
        totalCorrector
    end

    properties (Access = private)
        orientationVector
        dilatedOrientation
        mesh
        isCoherent
    end

    methods (Access = public)

        function obj = ConformalMappingComputer(cParams)
            obj.init(cParams);
            obj.isOrientationCoherent();
            obj.createInterpolator();
        end

        function phiV = compute(obj)
            obj.computeMappingWithSingularities();
            phiV = obj.phi;
        end

        function plot(obj)
            for iDim = 1:obj.mesh.ndim
                obj.phi{iDim}.plotContour();
            end
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.mesh               = cParams.mesh;
            obj.orientationVector  = cParams.orientationVector;
            obj.dilatedOrientation = cParams.dilatedOrientation;
        end

        function isOrientationCoherent(obj)
            s.mesh        = obj.mesh;
            s.orientation = obj.orientationVector{1}.project('P1D');
            c = CoherentOrientationSelector(s);
            isC = c.isOrientationCoherent();
            obj.isCoherent = isC;
        end        

        function createInterpolator(obj)
            s.mesh       = obj.mesh;
            s.isCoherent = obj.isCoherent;
            s = SymmetricContMapCondition(s);
            sC = s.computeCondition();
            obj.interpolator = sC;
        end

        function computeMappings(obj)
            nDim = obj.mesh.ndim;
            for iDim = 1:nDim
                bI    = obj.dilatedOrientation{iDim}.fValues;
                phiD  = obj.computeMapping(bI);
                phiV(iDim,:,:) = phiD.fValues;
            end
            s.fValues = phiV;
            s.mesh    = obj.mesh;
            phiF = P1DiscontinuousFunction(s);
            obj.phiMapping = phiF;
        end

        function computeTotalCorrector(obj)
           s.mesh = obj.mesh;
           s.isCoherent         = obj.isCoherent;
           s.dilatedOrientation = obj.dilatedOrientation;
           s.interpolator = obj.interpolator;
           s.phiMapping   = obj.phiMapping;
           tC = TotalCorrectorComputer(s);           
           obj.totalCorrector = tC.compute();
        end

        function computeMappingWithSingularities(obj)
            obj.computeMappings();
            obj.computeTotalCorrector();
            psiTs = P1DiscontinuousFunction.sum(obj.phiMapping,obj.totalCorrector);
            obj.phi = psiTs;
        end

        function phi = computeMapping(obj,fValues)
            s.fValues = fValues';
            s.mesh    = obj.mesh;
            s.rhsType = 'ShapeDerivative';
            s.interpolator = obj.interpolator;
            problem   = MinimumDiscGradFieldWithVectorInL2(s);
            phi = problem.solve();
        end

    end

end
