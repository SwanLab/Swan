classdef RHSintegrator < handle

    properties (Access = private)
        fGauss
        xGauss
        mesh
        type
        quadOrder
        
        quadrature
        nunknPerField
    end

    methods (Access = public, Static)
        
        function obj = create(s)
            f = RHSintegratorFactory();
            obj = f.create(s);
        end

    end
    
    methods (Access = public)
        
        function obj = RHSintegrator(cParams)
            obj.init(cParams);
            obj.computeQuadrature();
        end
        
        function int = integrate(obj)
            fG     = obj.fGauss;
            dV     = obj.computeDvolume();
%             fdV    = (fG.*dV);
            shapes = obj.computeShapeFunctions();
            nnode  = size(shapes,1);
            nelem  = size(shapes,2);
            int = zeros(nnode,nelem);
            for igaus = 1:obj.quadrature.ngaus
                nunkn = obj.nunknPerField;
                for iField = 1:nunkn % Apparently it really is always 1
                    fdv = fG(igaus,:,iField).*dV(igaus,:);
                    shape = shapes(:, :, igaus);
                    int = int + bsxfun(@times,shape,fdv);
                end
            end
            int = transpose(int);
        end
        
        function int = integrateWithShapeDerivative(obj)
            fG      = obj.fGauss;
            dV      = obj.computeDvolume();
            grad    = obj.computeGrad();
            nnode   = size(grad,2);
            ndim    = size(grad,1);
            nelem   = obj.mesh.nelem;
            int = zeros(nnode,nelem);
            for igaus = 1:obj.quadrature.ngaus
                for idime = 1:ndim
                    fI     = squeezeParticular(fG(idime,igaus,:),1);
                    fdV    = (fI.*dV(igaus));
                    dShape = squeeze(grad(idime,:,:,igaus));
                    intI = bsxfun(@times,dShape,fdV);
                    int = int + intI;
                end
            end
            int = transpose(int);
        end

    end
    
    methods (Access = private)

        function grad = computeGrad(obj)
            m.type = obj.mesh.type;
            int = Interpolation.create(m,'LINEAR');
            int.computeShapeDeriv(obj.xGauss);
            s.mesh = obj.mesh;
            g = Geometry.create(s);
            g.computeGeometry(obj.quadrature,int);
            grad = g.dNdx;
        end

        function init(obj,cParams)
            obj.fGauss    = cParams.fGauss;
            obj.xGauss    = cParams.xGauss;
            obj.mesh      = cParams.mesh;
            obj.type      = cParams.type;
            obj.quadOrder = cParams.quadOrder;
            obj.nunknPerField = 1;
        end
        
        function q = computeQuadrature(obj)
           q = Quadrature.set(obj.mesh.type);
           q.computeQuadrature(obj.quadOrder);
           obj.quadrature = q;
        end
        
        function dV = computeDvolume(obj)
            q = obj.quadrature;
            dV = obj.mesh.computeDvolume(q);
        end
        
        function shapes = computeShapeFunctions(obj)
            m.type = obj.type;
            int = Interpolation.create(m,'LINEAR');
            int.computeShapeDeriv(obj.xGauss);
            shapes = permute(int.shape,[1 3 2]);
        end
        
    end
    
end