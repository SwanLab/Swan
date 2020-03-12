classdef HomogenizedVarComputerFromInterpolation ...
        < HomogenizedVarComputer
    
    properties (Access = private)
        interpolation
        material
        pdim
        nElem
    end
    
    methods (Access = public)
        
        function obj = HomogenizedVarComputerFromInterpolation(cParams)
            obj.init(cParams);
            obj.createMaterialInterpolation(cParams);
            obj.createMaterial();
            obj.designVariable = cParams.designVariable;
        end
        
        function computeCtensor(obj,rho)
            nGaus  = size(rho,2);  
            nStres = obj.material.nstre;
            obj.C  = zeros(nStres,nStres,obj.nElem,nGaus);
            obj.dC = zeros(nStres,nStres,nGaus,obj.nElem);
            for igaus = 1:nGaus           
                rhoV = rho(:,igaus);
                p = obj.interpolation.computeMatProp(rhoV);
                obj.C(:,:,:,igaus)  = obj.computeC(p.mu,p.kappa);
                obj.dC(:,:,igaus,:) = obj.computeC(p.dmu,p.dkappa);
            end
        end
        
        function computeDensity(obj,rho)
            obj.rho = rho;
            obj.drho = ones(size(rho));
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.nElem = cParams.nelem;
            obj.pdim  = cParams.dim;
        end
        
        function createMaterial(obj)
            s.pdim  = obj.pdim;            
            s.ptype = 'ELASTIC';
            obj.material = Material.create(s);            
        end
        
        function createMaterialInterpolation(obj,cParams)
            s = cParams.interpolationSettings;
            s.nElem = obj.nElem;
            s.dim  = obj.pdim;
            int = MaterialInterpolation.create(s);
            obj.interpolation = int;
        end
        
        function C = computeC(obj,mu,kappa)
            s.kappa = kappa;
            s.mu    = mu;
            obj.material.compute(s);
            C  = obj.material.C;            
        end
        
    end
    
end