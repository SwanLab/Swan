classdef P1DiscontinuousFunction < FeFunction
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        connec
        fNodes
        interpolation
        type
    end
    
    properties (Access = private)
    end
    
    methods (Access = public)
        
        function obj = P1DiscontinuousFunction(cParams)
            obj.init(cParams)
            obj.createInterpolation();
        end

        function fxV = evaluate(obj, xV)
            % Goal: develop this function
            func = obj.fValues;
            obj.interpolation.computeShapeDeriv(xV);
            shapes = obj.interpolation.shape;
            nNode  = size(shapes,1);
            nGaus  = size(shapes,2);
            nF     = size(func,1);
            nElem  = size(func,3);
            fxV = zeros(nF,nGaus,nElem);
            for kNode = 1:nNode
                shapeKJ = shapes(kNode,:,:);
                fKJ     = func(:,kNode,:);
                f = bsxfun(@times,shapeKJ,fKJ);
                fxV = fxV + f;
            end
        end

        function plot(obj, m)
            dim = 1;
            ndims   = size(obj.fValues, 1);
            nelem   = size(obj.connec, 1);
            nnodeEl = size(obj.connec, 2);
            mD = m.createDiscontinuousMesh();
            fD = reshape(obj.fValues, [ndims, nelem*nnodeEl])';
            x = mD.coord(:,1);
            y = mD.coord(:,2);
            z = fD(:,dim);
            figure()
            a = trisurf(mD.connec, x, y, z);
%             colorbar
            shading interp
            a.EdgeColor = [0 0 0];
            view(0,90)
        end

    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.fValues = cParams.fValues;
            obj.connec  = cParams.connec;
            obj.type    = cParams.type;
            obj.ndimf   = size(cParams.fValues,1);
        end

        function createInterpolation(obj)
            m.type = obj.type;
            obj.interpolation = Interpolation.create(m,'LINEAR');
        end
        
    end
    
end