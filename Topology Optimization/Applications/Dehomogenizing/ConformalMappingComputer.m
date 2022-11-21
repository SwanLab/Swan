classdef ConformalMappingComputer < handle
    
    properties (Access = public)
       phi
       dilation
       corrector
    end
    
    properties (Access = private)
       
    end
    
    properties (Access = private)
       orientation 
       mesh
    end
    
    methods (Access = public)
        
        function obj = ConformalMappingComputer(cParams)
            obj.init(cParams);
        end
        
        function phiV = compute(obj)
            obj.computeDilation();
            obj.computeMapping();
            obj.computeOrthogonalCorrector();
            phiV = obj.phi;
        end
        
        function plot(obj)
            obj.plotDilation();
            obj.plotMapping();
        end
        
    end
    
    methods (Access = private)
        
        function plotDilation(obj)
           obj.plotField(obj.dilation);
        end
        
        function plotMapping(obj)
           phi1 = obj.phi(:,1);
           phi2 = obj.phi(:,2);
           obj.plotContour((phi1)); 
           obj.plotContour((phi2));
        end

        function plotField(obj,z)
            figure()
            s.mesh  = obj.mesh;
            s.field = z;
            n = NodalFieldPlotter(s);
            n.plot();
            shading interp
        end
        
        function plotContour(obj,z)
            figure()
            x = obj.mesh.coord(:,1);
            y = obj.mesh.coord(:,2);
            [~,h] = tricontour(obj.mesh.connec,x,y,z,50);
            set(h,'LineWidth',5);
            colorbar
        end
        
        function init(obj,cParams)
            obj.orientation = cParams.theta;
            obj.mesh        = cParams.mesh;
        end
        
        function computeDilation(obj)
            s.theta = obj.orientation;
            s.mesh  = obj.mesh;
            d = DilationFieldComputer(s);
            obj.dilation = d.compute();
        end
        
        function computeMapping(obj)
           nDim = 2;
           nnod = obj.mesh.nnodes;
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
            s.fGauss  = obj.computeVectorComponentP0(idim);
            s.fValues = obj.computeVector(idim);
            s.mesh   = obj.mesh;
            s.rhsType = 'ShapeFunction';
            varProb  = MinimumDiscGradFieldWithVectorInL2(s);
           % varProb  = MinimumGradFieldWithVectorInL2(s);            
            phi = varProb.solve();
        end
        
        function b = computeVector(obj,idim)
           er = exp(obj.dilation);
           erCos = er.*cos(obj.orientation);
           erSin = er.*sin(obj.orientation);
           Q(1,1,:) = erCos;
           Q(1,2,:) = -erSin;
           Q(2,1,:) = erSin;
           Q(2,2,:) = erCos;
           b = squeezeParticular(Q(:,idim,:),2);
        end
        
        function fG = computeVectorComponentP0(obj,idim)
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
            %%%%% HEREEEE!!!!! Integrate with more gauss points b
        end
        
        function computeOrthogonalCorrector(obj)
            s.mesh               = obj.mesh;
            s.orientation        = obj.orientation;
            o = OrthogonalCorrectorComputer(s);
            o.compute();
            o.plot();                        
        end
    
    end
    
end