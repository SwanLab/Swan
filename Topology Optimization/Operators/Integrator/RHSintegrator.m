classdef RHSintegrator < handle

    properties (Access = private)
        fGauss
        xGauss        
        mesh
        type
        quadOrder
        
        quadrature        
    end
    
    methods (Access = public)
        
        function obj = RHSintegrator(cParams)
            obj.init(cParams)   
            obj.computeQuadrature();
        end
        
        function int = integrate(obj)
            fG     = obj.fGauss;
            dV     = obj.computeDvolume();
            fdV    = (fG.*dV);
            shapes = obj.computeShapeFunctions();
            nnode  = size(shapes,1);
            nelem  = size(shapes,2);
            int = zeros(nnode,nelem);
            for igaus = 1:obj.quadrature.ngaus
                fdv   = fdV(igaus,:);
                shape = shapes(:,:,igaus);
                int = int + bsxfun(@times,shape,fdv);
            end
            int = transpose(int);            
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.fGauss    = cParams.fGauss;
            obj.xGauss    = cParams.xGauss;
            obj.mesh      = cParams.mesh;
            obj.type      = cParams.type;
            obj.quadOrder = cParams.quadOrder;
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