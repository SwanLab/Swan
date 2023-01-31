classdef OrthogonalCorrectorComputer < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        
        shiftingValue
        orthogonalCorrectorValue
    end
    
    properties (Access = private)
        mesh
        interpolator        
        correctorValue               
    end
    
    methods (Access = public)
        
        function obj = OrthogonalCorrectorComputer(cParams)
            obj.init(cParams)            
        end
        
        function c = compute(obj)
            obj.createShifting();
            obj.createOrthogonalCorrector();
            c = obj.orthogonalCorrectorValue;
        end
        
        function plot(obj)
            obj.plotFieldDG((obj.orthogonalCorrectorValue))
            obj.plotFieldDG((obj.shiftingValue))
            
            
            figure()
            m = obj.mesh.createDiscontinousMesh();
            x = m.coord(:,1);
            y = m.coord(:,2);
            z = abs(obj.shiftingValue');
            %figure()
            tricontour(m.connec,x,y,z,linspace(min(z(:)),max(z(:)),30))
         %   view(0,90)
        %    colorbar
        %    shading interp
        end
        
        function plotFieldDG(obj,f)
            figure()
            s.mesh  = obj.mesh.createDiscontinousMesh();
            s.field = transpose(f);
            n = NodalFieldPlotter(s);
            n.plot();
            shading interp             
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.mesh               = cParams.mesh;
            obj.interpolator       = cParams.interpolator;
            obj.correctorValue     = cParams.correctorValue;
        end

        function createShifting(obj)
            s.mesh     = obj.mesh;
            s.fGauss   = obj.createDiscontinousGradient(obj.correctorValue);
            s.rhsType = 'ShapeDerivative';
            s.interpolator = obj.interpolator;
            m = MinimumDiscGradFieldWithVectorInL2(s);
            f = m.solve();
            obj.shiftingValue = f;
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
        
       function createOrthogonalCorrector(obj)
            phi = obj.correctorValue;
            fD  = obj.shiftingValue;
            phi = phi -fD; 
            obj.orthogonalCorrectorValue = phi;
       end
       

    end
    
end