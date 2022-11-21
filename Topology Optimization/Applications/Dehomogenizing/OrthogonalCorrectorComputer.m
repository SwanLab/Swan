classdef OrthogonalCorrectorComputer < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        correctorValue
        shiftingValue
        orthogonalCorrectorValue
    end
    
    properties (Access = private)
        mesh
        orientation        
    end
    
    methods (Access = public)
        
        function obj = OrthogonalCorrectorComputer(cParams)
            obj.init(cParams)            
        end
        
        function compute(obj)
            obj.createCorrector();
            obj.createShifting();
            obj.createOrthogonalCorrector();
        end
        
        function plot(obj)
            phi = obj.orthogonalCorrectorValue;
            figure()
            s.mesh  = obj.mesh.createDiscontinousMesh();
            s.field = transpose(phi);
            n = NodalFieldPlotter(s);
            n.plot();
            shading interp             
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.mesh               = cParams.mesh;
            obj.orientation        = cParams.orientation;
        end
        
        function createCorrector(obj)   
            s.mesh               = obj.mesh;            
            s.orientation        = obj.orientation;
            c = CorrectorComputer(s);
            cV = c.compute();
            c.plot()                        
            obj.correctorValue = cV;
        end
        
        function createShifting(obj)
            s.mesh     = obj.mesh;
            s.fGauss   = obj.createCorrectorGradient();
            s.fValues  = obj.orientation';                        
            s.rhsType = 'ShapeDerivative';
            m = MinimumDiscGradFieldWithVectorInL2(s);
            f = m.solve();
            fD = obj.createDiscontinousField(f);
            fD  = reshape(fD',3,[])';  %%%Ehhhhhh!  
            obj.shiftingValue = fD;
        end
        
        function fGauss = createCorrectorGradient(obj) %%%Ehhhhh
            cV = obj.correctorValue;
            q = Quadrature.set(obj.mesh.type);
            q.computeQuadrature('LINEAR');              
            int = Interpolation.create(obj.mesh,obj.mesh.interpolation.order);
            int.computeShapeDeriv(q.posgp); 
            dN = int.deriv;
            nDim = size(dN,2);
            fGauss = zeros(nDim,q.ngaus,obj.mesh.nelem);
            for igaus = 1:q.ngaus
                for idim = 1:nDim
                    grad = dN(igaus,idim)*cV(:,idim);
                    fG(:,1) = squeeze(fGauss(idim,igaus,:));
                    fGauss(idim,igaus,:) = fG + grad;
                end
            end            
        end        
        
       function fD = createDiscontinousField(obj,fValues)
            s.connec = obj.mesh.connec;
            s.type   = obj.mesh.type;
            s.fNodes = fValues;
            f = FeFunction(s);            
            fD = f.computeDiscontinousField();
       end     
       

       function createOrthogonalCorrector(obj)
            phi = obj.correctorValue;
            fD  = obj.shiftingValue;
            phi = phi -fD; 
            obj.orthogonalCorrectorValue = phi;
       end
        
    end
    
end