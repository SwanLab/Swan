classdef HomogenizedVarComputerFromInterpolation ...
        < HomogenizedVarComputer
    
    properties (Access = private)
        interpolation
        material
    end
    
    methods (Access = public)
        
        function obj = HomogenizedVarComputerFromInterpolation(cParams)
            obj.init(cParams);
            obj.createMaterialInterpolation(cParams);
            obj.createMaterial(cParams);
        end
        
        function computeCtensor(obj,x)
            rho = x{1};
            nvar = 1;
            nGaus  = size(rho,2);
            nElem  = size(rho,1);
            nStres = obj.material.nstre;
            obj.C  = zeros(nStres,nStres,nElem,nGaus);
            obj.dC = zeros(nStres,nStres,nvar,nElem,nGaus);
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
            nElem  = size(rho,1);
            nStres = obj.material.nstre;
            obj.Pp   = zeros(nStres,nStres,nElem,nGaus);
            obj.dPp = zeros(nStres,nStres,nvar,nElem,nGaus);
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
        end
        
        function createMaterial(obj,cParams)
            if isprop(cParams,'material')
                s.pdim  = cParams.dim;
                s.ptype = 'ELASTIC';
              obj.material = Material.create(s);
            else
              obj.material = cParams.material;
            end
        end
        
        function createMaterialInterpolation(obj,cParams)
            s = cParams.interpolationSettings;
            if isprop(s,'interpolation')
                s.nElem = obj.mesh.nElem;
                s.dim   = cParams.dim;
                int = MaterialInterpolation.create(s);
                obj.interpolation = int;
            else
                obj.interpolation = s.interpolation;
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