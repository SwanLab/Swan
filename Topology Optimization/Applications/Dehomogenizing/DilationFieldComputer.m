classdef DilationFieldComputer < handle
    
    properties (Access = private)
       dilation 
    end
    
    properties (Access = private)
       theta 
       mesh
    end
    
    methods (Access = public)
        
        function obj = DilationFieldComputer(cParams)
            obj.init(cParams)
        end
        
        function d = compute(obj)
            obj.computeDilationField();
            d = obj.dilation; 
        end
        
        function plot(obj)
            figure()
            s.mesh  = obj.mesh;
            s.field = obj.dilation;
            n = NodalFieldPlotter(s);
            n.plot();
            shading interp            
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.theta = cParams.theta;
            obj.mesh  = cParams.mesh;
        end
               
        function computeDilationField(obj)
            s.fGauss = obj.computeThetaGradient();
            s.mesh   = obj.mesh;
            varProb  = MinimumGradFieldWithVectorInL2(s);
            r = varProb.solve();
            obj.dilation = r;
        end
        
        function gradT = computeThetaGradient(obj)
            q = Quadrature.set(obj.mesh.type);
            q.computeQuadrature('CUBIC');
            m.type = obj.mesh.type;
            int = Interpolation.create(m,'LINEAR');
            int.computeShapeDeriv(q.posgp);
            s.mesh = obj.mesh;
            g = Geometry.create(s);
            g.computeGeometry(q,int);
          %  grad = g.dNdx;
            nodes = obj.mesh.connec;


            alpha = obj.theta;
            a1(:,1) = cos(alpha);
            a1(:,2) = sin(alpha);

            a2(:,1) = -sin(alpha);
            a2(:,2) = cos(alpha);

            dN = g.dNdx;

            s.connec  = obj.mesh.connec;
            s.type    = obj.mesh.type;
            s.fValues = a1;            
            a1F = P1Function(s);
            a1FV = a1F.evaluate(q.posgp); 


            s.connec  = obj.mesh.connec;
            s.type    = obj.mesh.type;
            s.fValues = a2;            
            a2F = P1Function(s);
            a2FV = a2F.evaluate(q.posgp);            

            da1 = obj.computeDivergence(a1F,q,dN);
            da2 = obj.computeDivergence(a2F,q,dN);            
            
            da1a1 = bsxfun(@times,da1,a1FV);
            da2a2 = bsxfun(@times,da2,a2FV);
          
            
            gradT3 = -da2a2 - da1a1;
         %   norm(gradT3(:) - gradT2(:))/norm(gradT2(:))

            gradT = gradT3;
        end
        
        function da = computeDivergence(obj,f,q,dN)
            fV = f.fValues;
            nodes = obj.mesh.connec;
            nNode = obj.mesh.nnodeElem;
            nDim  = obj.mesh.ndim;
            da = zeros(1,q.ngaus,obj.mesh.nelem);
            for igaus = 1:q.ngaus
                for kNode = 1:nNode
                    nodeK = nodes(:,kNode);
                    for rDim = 1:nDim
                        dNkr = squeeze(dN(rDim,kNode,:,igaus));
                        fkr = fV(nodeK,rDim);
                        int(1,1,:) = dNkr.*fkr;
                        da(1,igaus,:) = da(1,igaus,:) + int(1,1,:);
                    end
                end
            end
       end
        



    end

    
end