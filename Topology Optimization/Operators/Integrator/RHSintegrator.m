classdef RHSintegrator < handle

    properties (Access = private)
        fNodal
        backgroundMesh
        mesh
        feMesh
        xGauss
        quadrature
    end
    
    methods (Access = public)
        
        function obj = RHSintegrator(cParams)
            obj.init(cParams)            
        end
        
        function int = integrate(obj)
            Fgauss  = obj.computeFinGaussPoints();
            dV      = obj.computeDvolume();
            fdV     = (Fgauss.*dV);
            shapes  = obj.computeShapeFunctions();
            int = obj.initIntegrand();
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
            obj.fNodal         = cParams.fNodal;
            obj.backgroundMesh = cParams.backgroundMesh;
            obj.xGauss         = cParams.xGauss;
            obj.mesh           = cParams.mesh;
            obj.quadrature     = cParams.quadrature;
            obj.feMesh = cParams.feMesh;
        end
        
        function fG = computeFinGaussPoints(obj)
            f  = obj.createFeFunction();
            fG = f.interpolateFunction(obj.xGauss);
            fG = permute(fG,[2 3 1]);            
        end        
        
        function f = createFeFunction(obj)            
            m = obj.feMesh;
            s.fNodes = obj.fNodal;
            s.mesh = m;
            f = FeFunction(s);
        end             
        
        function dV = computeDvolume(obj)
            q  = obj.quadrature;
            dV = obj.mesh.computeDvolume(q);            
        end        
        
        function shapes = computeShapeFunctions(obj)
            m = obj.backgroundMesh;
            int = Interpolation.create(m,'LINEAR');
            int.computeShapeDeriv(obj.xGauss);
            shapes = permute(int.shape,[1 3 2]);            
        end          
        
        function int = initIntegrand(obj)
            nelem = obj.mesh.nelem;
            nnode = obj.backgroundMesh.nnode;
            int = zeros(nnode,nelem);            
        end          
        
    end
    
end