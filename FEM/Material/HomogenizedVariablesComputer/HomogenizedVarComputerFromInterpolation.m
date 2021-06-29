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
        end
        
        function computeCtensor(obj,x)
            rho = x{1};
            nvar = 1;
            nGaus  = size(rho,2);  
            nStres = obj.material.nstre;
            obj.C  = zeros(nStres,nStres,obj.nElem,nGaus);
            obj.dC = zeros(nStres,nStres,nvar,obj.nElem,nGaus);
            for igaus = 1:nGaus           
                rhoV = rho(:,igaus);
                p = obj.interpolation.computeMatProp(rhoV);
                obj.C(:,:,:,igaus)  = obj.computeC(p.mu,p.kappa);
                for ivar = 1:nvar
                    obj.dC(:,:,1,:,igaus) = obj.computeC(p.dmu,p.dkappa);
                end
            end
            obj.Cref = obj.C;
            obj.dCref = obj.dC;
        end
        
        function computeDensity(obj,x)
            rho = x{1};            
            obj.rho = rho;
            obj.drho{1} = ones(size(rho));
        end
        
        function computePtensor(obj,x,pNorm)
            rho = x{1};
            nvar = 1;
            nGaus  = size(rho,2);  
            nStres = obj.material.nstre;
            obj.Pp   = zeros(nStres,nStres,obj.nElem,nGaus);
            obj.dPp = zeros(nStres,nStres,nvar,obj.nElem,nGaus);
            q = 2.5;
            Pr  = @(rho) 1./rho.^(q);
            dPr = @(rho) (-q)*(1./(rho.^(q+1)));
            for iStres = 1:nStres
                obj.Pp(iStres,iStres,:,:)  = Pr(rho);
                for ivar = 1:nvar
                    obj.dPp(iStres,iStres,ivar,:,:) = dPr(rho)';
                end
            end            
        end
        
        function fP = addPrintableVariables(obj,x)            
            fP{1}.value = obj.rho;
        end
        
        function fP = createPrintVariables(obj)
            fP{1}.type = 'ScalarGauss';            
            fP{1}.name = 'DensityGauss';            
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
            if ~isempty(s.typeOfMaterial)
                s.nElem = obj.nElem;
                s.dim  = obj.pdim;
                int = MaterialInterpolation.create(s);
                obj.interpolation = int;
            end
        end
        
        function C = computeC(obj,mu,kappa)
            s.kappa = kappa;
            s.mu    = mu;
            obj.material.compute(s);
            C  = obj.material.C;            
        end
        
    end
    
end