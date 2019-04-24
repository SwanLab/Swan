classdef HomogenizedVarComputerFromVademecum ...
        < HomogenizedVarComputer
    
    properties (Access = public)
        Ctensor
        density
    end
    
    
    methods (Access = public)
        
        function obj = HomogenizedVarComputerFromVademecum(cParams)
            s.fileName = cParams.fileName;
            v = VademecumVariablesLoader(s);
            v.load();
            obj.Ctensor = v.Ctensor;
            obj.density = v.density;
        end
        
        function computeCtensor(obj,rho)
            rho = max(0.001,min(rho,0.999));
            mx = max(0.01,min(sqrt(1-rho),0.99));
            my = max(0.01,min(sqrt(1-rho),0.99));            
            [c,dc] = obj.Ctensor.compute([mx,my]);
            obj.dC = zeros(3,3,size(dc,3));
            obj.C = c;
            for i = 1:3
                for j = 1:3
                    obj.dC(i,j,:) = squeeze(-dc(i,j,:,1))./my;
                end
            end                                    
        end
        
        function computeDensity(obj,rho)
            rho = max(0.001,min(rho,0.999));
            mx = max(0.01,min(sqrt(1-rho),0.99));
            my = max(0.01,min(sqrt(1-rho),0.99));            
            [rho,drho] = obj.density.compute([mx,my]);
            obj.rho = rho;
            obj.drho = ones(size(rho));
        end
        
    end
    
end