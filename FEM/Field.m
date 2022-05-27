classdef Field < handle

    properties (Access = public)
        dim
        connec
        coord % should be removed
        geometry % moved to Mesh
        fefnc % merge
        interpolation
        xGauss % no
        quadrature % perhaps private
    end

    properties (Access = private)
        mesh
        ndimf
        scale
%         quadrature
%         interpolation
        interpTranslator % new class to do that from Field
        interpolationOrder
        quadratureOrder
    end

    methods (Access = public)

        function obj = Field(cParams)
            obj.init(cParams)
            obj.createQuadrature();
            obj.createInterpolation();
            obj.createGeometry();
            obj.computeXGauss();
            obj.updateInputMismatch();
            obj.computeDimensions();
        end

        function newBC = translateBoundaryConditions(obj, inBC)
            newBC = obj.interpTranslator.updateBoundaryConditions(inBC);
            
        end

        function newBC = adaptBoundaryConditions(obj, inBC)
            fn = fieldnames(inBC);
            for k=1:numel(fn)
                param = inBC.(fn{k});
                if isstruct(param)
                    switch param.domain
                        case 'Border'
                        x = obj.coord(:,1);
                        y = obj.coord(:,2);
                        idx = find(x == 0 | x == 1 | y == 0 | y == 1);
                        one = ones(length(idx),1);
                        newDirich = [idx, one, one*param.value;
                                     idx, 2*one, one*param.value];
                        newBC.dirichlet = sortrows(newDirich);
                    end
                end
            end
        end

    end

    methods (Access = private)

        function init(obj, cParams)
            obj.mesh               = cParams.mesh;
            obj.ndimf              = cParams.ndimf;
            obj.scale              = cParams.scale;
            obj.interpolationOrder = cParams.interpolationOrder;
            if isfield(cParams, 'quadratureOrder')
                obj.quadratureOrder = cParams.quadratureOrder;
            else
                obj.quadratureOrder = cParams.interpolationOrder;
            end

        end

        function createQuadrature(obj)
            quad = Quadrature.set(obj.mesh.type);
            quad.computeQuadrature(obj.quadratureOrder);
            obj.quadrature = quad;
        end

        function createInterpolation(obj)
            int = Interpolation.create(obj.mesh,obj.interpolationOrder);
            int.computeShapeDeriv(obj.quadrature.posgp);
            obj.interpolation = int;
        end
       
        function createGeometry(obj)
            q = obj.quadrature;
            int = obj.interpolation;
            int.computeShapeDeriv(q.posgp);
            s.mesh = obj.mesh;
            g = Geometry.create(s);
            g.computeGeometry(q,int);
            obj.geometry = g;
        end

        function computeXGauss(obj)
            xV = obj.quadrature.posgp;
            xG = obj.mesh.computeXgauss(xV);
            obj.xGauss = xG;
        end

        function updateInputMismatch(obj) % perhaps should be renamed
%             if ~isequal(obj.interpolation.order, 'LINEAR')
                s.mesh          = obj.mesh;
                s.interpolation = obj.interpolation;
                obj.interpTranslator = InterpolationTranslator(s);
                obj.coord  = obj.interpTranslator.coord;
                obj.connec = obj.interpTranslator.connec;
%             end
        end

        function computeDimensions(obj)
            s.coord         = obj.coord;
            s.ndimf         = obj.ndimf;
            s.interpolation = obj.interpolation;
            d = FieldDimensions(s);
            d.compute();
            obj.dim = d;
        end

    end

end