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

        function computeMappings(obj)
            s.mesh              = obj.mesh;
            s.orientationVector  = obj.orientationVector;
            s.dilatedOrientation = obj.dilatedOrientation;           
            mC   = MappingComputer(s);
            phiM = mC.compute();
            obj.phiMapping = phiM;
        end

        function computeTotalCorrector(obj)
           s.mesh = obj.mesh;
           s.orientationVector  = obj.orientationVector;
           s.dilatedOrientation = obj.dilatedOrientation;
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

    end

end
