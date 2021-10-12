classdef ConformalMappingComputer < handle
    
    properties (Access = public)
       phi
       dilation
    end
    
    properties (Access = private)
       
    end
    
    properties (Access = private)
       theta 
       mesh        
    end
    
    methods (Access = public)
        
        function obj = ConformalMappingComputer(cParams)
            obj.init(cParams);            
        end
        
        function phiV = compute(obj)
            obj.computeDilation();
            obj.computeMapping();
            phiV = obj.phi;
        end
        
        function plot(obj)
            obj.plotDilation();
            obj.plotMapping();
        end
        
    end
    
    methods (Access = private)
        
        function plotDilation(obj)
            figure()
            x = obj.mesh.coord(:,1);
            y = obj.mesh.coord(:,2);
            z = obj.dilation;
            F = scatteredInterpolant(x,y,z);
            n = 100;
            [xv,yv] = meshgrid(min(x):1/n:max(x),min(y):1/n:max(y));
            zv = F(xv,yv);
            %zv = griddata(x,y,z,xv,yv);
            surf(xv,yv,zv);
            view(0,90)
            colorbar
            shading interp
        end
        
        function plotMapping(obj)
           phi1 = obj.phi(:,1);
           phi2 = obj.phi(:,2);
           obj.plotContour(phi1); 
           obj.plotContour(phi2);
        end
        
        function plotContour(obj,z)
            figure()
            x = obj.mesh.coord(:,1);
            y = obj.mesh.coord(:,2);
            F = scatteredInterpolant(x,y,z);
            n  =100;
            [xv,yv] = meshgrid(min(x):1/n:max(x),min(y):1/n:max(y));
            zv = F(xv,yv);
            %zv = griddata(x,y,z,xv,yv);
            contour(xv,yv,zv,50,'LineWidth',5);
            colorbar
        end
        
        function init(obj,cParams)
            obj.theta    = cParams.theta;
            obj.mesh     = cParams.mesh;
        end
        
        function computeDilation(obj)
            s.theta = obj.theta;
            s.mesh  = obj.mesh;
            d = DilationFieldComputer(s);
            obj.dilation = d.compute();           
        end
        
        function computeMapping(obj)
           nDim = 2;
           nnod = obj.mesh.npnod;
           phiV = zeros(nnod,nDim);
           for iDim = 1:nDim
             phiV(:,iDim) = obj.computeComponent(iDim);
           end
           obj.phi = phiV;
        end

        function computeSecondComponent(obj)
           iDim = 2; 
           obj.phi(:,iDim) = obj.computeComponent(iDim);
        end
              
        function phi = computeComponent(obj,idim)
            s.fGauss = obj.computeVectorComponent(idim);
            s.mesh   = obj.mesh;
            varProb  = MinimumGradFieldWithVectorInL2(s);            
            phi = varProb.solve();
        end
        
        function b = computeVector(obj,idim)
           er = exp(obj.dilation); 
           erCos = er.*cos(obj.theta);
           erSin = er.*sin(obj.theta);
           Q(1,1,:) = erCos;
           Q(1,2,:) = -erSin;
           Q(2,1,:) = erSin;
           Q(2,2,:) = erCos;
           b = squeezeParticular(Q(:,idim,:),2);            
        end
        
        function fG = computeVectorComponent(obj,idim)
            b = obj.computeVector(idim);            
            q = Quadrature.set(obj.mesh.type);
            q.computeQuadrature('LINEAR');
            xGauss = q.posgp;
            fG = zeros(obj.mesh.ndim,q.ngaus,obj.mesh.nelem);
            for idim = 1:obj.mesh.ndim
                s.fNodes = b(idim,:)';
                s.connec = obj.mesh.connec;
                s.type   = obj.mesh.type;
                f = FeFunction(s);
                fG(idim,:,:) = f.interpolateFunction(xGauss);            
            end     
        end        
        
    end
    
end