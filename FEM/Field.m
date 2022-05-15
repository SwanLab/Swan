classdef Field < handle

    properties (Access = public)
        dim
        connec
        coord
        geometry
        boundaryConditions
        inputBC %private
        fefnc % 
        interpolation
    end

    properties (Access = private)
        mesh
        ndimf
        scale
        quadrature
%         interpolation
        interpTranslator
        interpolationOrder
    end

    methods (Access = public)

        function obj = Field(cParams)
            obj.init(cParams)
            obj.createQuadrature();
            obj.createInterpolation();
            obj.createGeometry();
            obj.updateInputMismatch();
            obj.computeDimensions();
            obj.translateBoundaryConditions();
            obj.createBoundaryConditions();
        end

    end

    methods (Access = private)

        function init(obj, cParams)
            obj.mesh               = cParams.mesh;
            obj.ndimf              = cParams.ndimf;
            obj.scale              = cParams.scale;
            obj.inputBC            = cParams.inputBC;
            obj.interpolationOrder = cParams.interpolationOrder;
        end

        function createQuadrature(obj)
            quad = Quadrature.set(obj.mesh.type);
            quad.computeQuadrature(obj.interpolationOrder);
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

        function updateInputMismatch(obj) % perhaps should be renamed
%             if ~isequal(obj.interpolation.order, 'LINEAR')
                s.mesh          = obj.mesh;
                s.interpolation = obj.interpolation;
                s.inputBC       = obj.inputBC;
                obj.interpTranslator = InterpolationTranslator(s);
                obj.coord  = obj.interpTranslator.coord;
                obj.connec = obj.interpTranslator.connec;
                obj.inputBC = obj.interpTranslator.inputBC;
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

        function translateBoundaryConditions(obj)
            inBC = obj.inputBC;
            fn = fieldnames(inBC);
            for k=1:numel(fn)
                param = inBC.(fn{k});
                if isstruct(param)
                    disp('we got ourselves a struct')
                    if (param.domain == 'Border')
                        x = obj.coord(:,1);
                        y = obj.coord(:,2);
                        idx = find(x == 0 | x == 1 | y == 0 | y == 1);
                        one = ones(length(idx),1);
                        newDirich = [idx, one, one*param.value;
                                     idx, 2*one, one*param.value];
                        obj.inputBC.dirichlet = newDirich;
                    end
                end
            end
        end

        function createBoundaryConditions(obj)
            s.dim   = obj.dim;
            s.mesh  = obj.mesh;
            s.scale = obj.scale;
            s.bc    = obj.inputBC;
%             s.ndofs = obj.dim.ndofs;
            bc = BoundaryConditions(s);
            bc.compute();
            obj.boundaryConditions = bc;
        end

    end

end