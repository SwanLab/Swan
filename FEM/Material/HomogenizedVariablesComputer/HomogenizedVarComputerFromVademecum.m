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
        
        function computeCtensor(obj,rho,princDir)
            R = obj.computeRotatorMatrix(princDir);
            rho = max(0.001,min(rho,0.999));
            mx = max(0.01,min(sqrt(1-rho),0.99));
            my = max(0.01,min(sqrt(1-rho),0.99));
            [c,dc] = obj.Ctensor.compute([mx,my]);
            obj.dC = zeros(3,3,size(dc,3));
            for i = 1:3
                for j = 1:3
                    dct = squeeze(-dc(i,j,:,1))./my;
                    obj.dC(i,j,:) = dct;
                end
            end
 
            cr = zeros(size(c));
            dCr = zeros(size(obj.dC));
            for i = 1:3
                for j= 1:3
                    for m = 1:3
                        for n = 1:3
                            Cij  = R(m,i,:).*c(m,n,:).*R(n,j,:);
                            dCij = R(m,i,:).*obj.dC(m,n,:).*R(n,j,:);
                            cr(i,j,:)  = cr(i,j,:) + Cij;
                            dCr(i,j,:) = dCr(i,j,:) + dCij;
                        end
                    end                    
                end
            end
            obj.C = cr;
            obj.dC = dCr;
%            obj.C = c;
            
        end
        
        function R = computeRotatorMatrix(obj,dir)
            angle = squeeze(acos(dir(1,1,:)));
            Rs = obj.createSymbolicRotationMatrix();
            R = zeros(3,3,length(angle));
            for i = 1:3
                for j = 1:3
                    r = matlabFunction(Rs(i,j));
                    R(i,j,:) = r(angle);
                end
            end
        end
        
        function a = createSymbolicRotationMatrix(obj)
            S = StressPlaneStressVoigtTensor();            
            alpha = sym('alpha');
            v = Vector3D();
            v.setValue([0 0 1]);
            factory = RotatorFactory();
            rotator = factory.create(S,alpha,v);
            a = rotator.rotationMatrix();            
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