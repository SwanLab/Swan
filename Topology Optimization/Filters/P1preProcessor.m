classdef P1preProcessor < handle

    properties (Access = private)
        mesh
        quadratureOrder
    end

    properties (Access = public)
        quadrature
        interp
        geometry
    end

    methods (Access = public)
        function obj = P1preProcessor(cParams)
            obj.init(cParams);
        end

        function preProcess(obj)
            obj.createQuadrature();
            obj.createInterpolation();
            obj.createGeometry();
        end
    end

    methods (Access = private)

        function init(obj,cParams)
            obj.mesh = cParams.mesh;
            obj.quadratureOrder = cParams.quadratureOrder;
        end

        function createInterpolation(obj)
            obj.interp = Interpolation.create(obj.mesh,'LINEAR');
        end

        function createGeometry(obj)
            s.mesh = obj.mesh;
            obj.geometry = Geometry.create(s);
            obj.geometry.computeGeometry(obj.quadrature,obj.interp);
        end

        function createQuadrature(obj)
            obj.quadrature = Quadrature.set(obj.mesh.type);
            obj.quadrature.computeQuadrature(obj.quadratureOrder);
        end

    end

end