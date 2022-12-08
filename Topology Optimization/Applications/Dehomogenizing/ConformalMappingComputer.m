classdef ConformalMappingComputer < handle
    
    properties (Access = public)
       phi
       dilation
       interpolator
       singularityCoord
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
            obj.createInterpolator();
            obj.computeDilation();
            obj.computeSingularities();
            obj.computeMappingWithSingularities();
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
            m = obj.mesh.createDiscontinousMesh;
            x = m.coord(:,1);
            y = m.coord(:,2);
            [~,h] = tricontour(m.connec,x,y,z,50);
            set(h,'LineWidth',5);
            colorbar
        end
        
        function init(obj,cParams)
            obj.orientation = cParams.theta;
            obj.mesh        = cParams.mesh;
        end
        
        function createInterpolator(obj)
            s.meshCont    = obj.mesh;
            s.meshDisc    = obj.mesh.createDiscontinousMesh();
            s.orientation = obj.orientation;
            s = SymmetricContMapCondition(s);            
            sC = s.computeCondition();
            obj.interpolator = sC;            
        end        
        
        function computeDilation(obj)
            s.theta = obj.orientation;
            s.mesh  = obj.mesh;
            d = DilationFieldComputer(s);
            obj.dilation = d.compute();
        end
        
        
        function computeSingularities(obj)
            s.mesh        = obj.mesh;
            s.orientation = obj.computeOrientationVector(1)';
            sF = SingularitiesFinder(s);
            isS = sF.computeSingularElements();
            sF.plot();
            coordB = obj.mesh.computeBaricenter();
            coordB = transpose(coordB);
            sCoord =  coordB(isS,:);
            obj.singularityCoord = sCoord;
        end        
        
        function computeMappingWithSingularities(obj)
           nDim = 2;
           m = obj.mesh.createDiscontinousMesh;
           nnod = m.nnodes;
           phiV = zeros(nnod,nDim);
           for iDim = 1:nDim
              b = obj.computeOrientationVector(iDim);             
              bGauss = obj.computeOrientationVectorComponentP0(b);               
              phiI   = obj.computeMapping(bGauss);              
              if ~isempty(obj.singularityCoord)
                cr = obj.computeCorrector(b);
                oC = obj.computeOrthogonalCorrector(cr);
                coef = obj.computeCoeffs(oC,bGauss);
                psiT = zeros(size(oC));
                for iSing = length(coef)
                  psiT = psiT + coef(iSing)*oC(:,:,iSing);
                end 
                phiIc = phiI + psiT;
              else
                phiIc = phiI; 
              end              
              phiV(:,iDim) = reshape(phiIc',[],1);
           end
           obj.phi = abs(phiV);
        end
        
        function phi = computeMapping(obj,fGauss)
            s.fGauss  = fGauss;            
            s.mesh    = obj.mesh;
            s.rhsType = 'ShapeDerivative';
            s.interpolator = obj.interpolator;
            varProb   = MinimumDiscGradFieldWithVectorInL2(s);
           % varProb  = MinimumGradFieldWithVectorInL2(s);            
            phi = varProb.solve();
        end        
        
        function cV = computeCorrector(obj,b) 
            sCoord = obj.singularityCoord;
            s.mesh               = obj.mesh;            
            s.orientation        = b';
            s.singularityCoord = sCoord(end-1,:);
            c = CorrectorComputer(s);
            cV = c.compute();
            c.plot()                        
        end
        

        function b = computeOrientationVector(obj,idim)
           er = exp(obj.dilation);
           erCos = er.*cos(obj.orientation);
           erSin = er.*sin(obj.orientation);
           Q(1,1,:) = erCos;
           Q(1,2,:) = -erSin;
           Q(2,1,:) = erSin;
           Q(2,2,:) = erCos;
           b = squeezeParticular(Q(:,idim,:),2);
        end
        
        function bG = computeOrientationVectorComponentP0(obj,b)
            q = Quadrature.set(obj.mesh.type);
            q.computeQuadrature('QUADRATIC');
            xGauss = q.posgp;
            bG = zeros(obj.mesh.ndim,q.ngaus,obj.mesh.nelem);
            for idim = 1:obj.mesh.ndim
                s.fNodes = b(idim,:)';
                s.connec = obj.mesh.connec;
                s.type   = obj.mesh.type;
                f = FeFunction(s);
                bG(idim,:,:) = f.interpolateFunction(xGauss);
            end
            %%%%% HEREEEE!!!!! Integrate with more gauss points b
        end
        
        function c = computeOrthogonalCorrector(obj,c)
            s.mesh               = obj.mesh;
            s.correctorValue     = c;
            s.interpolator       = obj.interpolator;
            o = OrthogonalCorrectorComputer(s);
            c = o.compute();
           % o.plot();
        end
        
        function c = computeCoeffs(obj,psi,bGauss)
            In = obj.interpolator;             
            q = Quadrature.set(obj.mesh.type);
            q.computeQuadrature('QUADRATIC');   
            m = obj.mesh.createDiscontinousMesh();
            dV = obj.mesh.computeDvolume(q);
            nDim = m.ndim;
            nSing = size(psi,3);
            LHS = zeros(nSing,nSing);
            RHS = zeros(nSing,1);
            
            nnode = obj.mesh.nnodeElem;
            nElem = obj.mesh.nelem;
            dPsi = zeros(nDim,nnode,nElem,nSing);
            for iSing = 1:nSing
                dPsi(:,:,:,iSing) = obj.createDiscontinousGradient(psi(:,:,iSing)); 
            end
            
            
            for igaus = 1:q.ngaus
                dVG  = dV(igaus,:)';
                for idim = 1:nDim
                    bG = squeeze(bGauss(idim,igaus,:));
                    dPsiG = squeeze(dPsi(idim,igaus,:,:));
                    for iSing = 1:nSing
                        dPsiI = dPsiG(:,iSing);
                        for jSing = 1:nSing
                            dPsiJ = dPsiG(:,:,:,jSing);
                            lhs = sum(dPsiI.*dVG.*dPsiJ);
                            LHS(iSing,jSing) = LHS(iSing,jSing) + lhs;
                        end
                        rhs = sum(dPsiI.*dVG.*bG);
                        RHS(iSing) = RHS(iSing) + rhs;
                    end
                    
                end
                
            end
            c = LHS\RHS;
        end 
        
        function fGauss = createDiscontinousGradient(obj,field)%%%Ehhhhh            
            cV = field; 
            q = Quadrature.set(obj.mesh.type);
            q.computeQuadrature('QUADRATIC');              
            int = Interpolation.create(obj.mesh,obj.mesh.interpolation.order);
            int.computeShapeDeriv(q.posgp); 
            s.mesh = obj.mesh;
            g = Geometry.create(s);
            g.computeGeometry(q,int);
            dN = g.dNdx;            
            
            %dN = int.deriv;
            nDim = size(dN,1);
            nnode = size(dN,2);
            fGauss = zeros(nDim,q.ngaus,obj.mesh.nelem);
            for igaus = 1:q.ngaus
                for idim = 1:nDim
                    for iNode = 1:nnode
                    dNi = squeeze(dN(idim,iNode,:,igaus));
                    grad = dNi.*cV(:,iNode);
                    fG(:,1) = squeeze(fGauss(idim,igaus,:));
                    fGauss(idim,igaus,:) = fG + grad;
                    end
                end
            end            
        end           
                
    
    end
    
end