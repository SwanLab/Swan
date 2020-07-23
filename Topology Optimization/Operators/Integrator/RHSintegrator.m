classdef RHSintegrator < handle

    properties (Access = private)
        fNodal
        mesh
        type
        connec
        xGauss
        quadrature
    end
    
    methods (Access = public)
        
        function obj = RHSintegrator(cParams)
            obj.init(cParams)   
            obj.computeQuadrature();
        end
        
        function int = integrate(obj)
            Fgauss  = obj.computeFinGaussPoints();
            dV      = obj.computeDvolume();
            fdV     = (Fgauss.*dV);
            shapes  = obj.computeShapeFunctions();
            nnode   = size(shapes,1);
            nelem   = size(shapes,2);
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
            obj.fNodal = cParams.fNodal;
            obj.xGauss = cParams.xGauss;
            obj.mesh   = cParams.mesh;
            obj.type   = cParams.type;
            obj.connec = cParams.connec;
        end
        
        function q = computeQuadrature(obj)
           q = Quadrature.set(obj.mesh.type);
           q.computeQuadrature('LINEAR');
           obj.quadrature = q;
       end          
        
        function fG = computeFinGaussPoints(obj)
            f  = obj.createFeFunction();
            fG = f.interpolateFunction(obj.xGauss);
            fG = permute(fG,[2 3 1]);            
        end        
        
        function f = createFeFunction(obj)            
            s.fNodes = obj.fNodal;
            s.connec = obj.connec;
            s.type   = obj.type;
            f = FeFunction(s);
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