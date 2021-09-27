classdef DilationFieldComputer < handle
    
    properties (Access = private)
       theta 
       mesh
    end
    
    methods (Access = public)
        
        function obj = DilationFieldComputer(cParams)
            obj.init(cParams)            
        end
        
        function d = compute(obj)
            d = obj.computeDilationField();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.theta = cParams.theta;
            obj.mesh  = cParams.mesh;
        end
               
        function r = computeDilationField(obj)
            s.fGauss = obj.computeThetaGradient();
            s.mesh   = obj.mesh;
            varProb  = MinimumGradFieldWithVectorInL2(s);            
            r = varProb.solve();
        end
        
        function gradT = computeThetaGradient(obj)
            q = Quadrature.set(obj.mesh.type);
            q.computeQuadrature('LINEAR');
            m.type = obj.mesh.type;
            int = Interpolation.create(m,'LINEAR');
            int.computeShapeDeriv(q.posgp);
            s.mesh = obj.mesh;
            g = Geometry.create(s);
            g.computeGeometry(q,int);
            grad = g.cartd;
            nodes = obj.mesh.connec;
            f = obj.theta;
            gradF = zeros(obj.mesh.ndim,q.ngaus,obj.mesh.nelem);
            for igaus = q.ngaus
                for inode = 1:obj.mesh.nnode
                   nodeI = nodes(:,inode);
                   fI = f(nodeI);  
                   for idim = 1:obj.mesh.ndim
                    dN = squeeze(grad(idim,inode,:,igaus));    
                    gF = squeeze(gradF(idim,igaus,:)); 
                    gradF(idim,igaus,:) = gF + fI.*dN;%bsxfun(@times,dN,fI);
                   end
                end
            end
            gradT = zeros(size(gradF));
            gradT(1,:,:) = -gradF(2,:,:);
            gradT(2,:,:) = gradF(1,:,:);
        end        
        
    end
    
end