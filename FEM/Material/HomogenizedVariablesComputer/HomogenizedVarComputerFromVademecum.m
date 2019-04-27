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
        
        function computeCtensor(obj,x,princDir)
            R = obj.computeRotatorMatrix(princDir);
            %rho = max(0.001,min(rho,0.999));
            %mx = max(0.01,min(sqrt(1-rho),0.99));
            %my = max(0.01,min(sqrt(1-rho),0.99));
            %nv = length(x)/2;
            %mx = x(1:nv);
            %my = x(nv+1:2*nv,2);
            mx = x(:,1);
            my = x(:,2);
            [c,dc] = obj.Ctensor.compute([mx,my]);
            dc = permute(dc,[1 2 4 3]);
            % obj.dC = zeros(dc);
            %for i = 1:3
            %    for j = 1:3
            %        dct = squeeze(-dc(i,j,:,1))./my;
            %        obj.dC(i,j,:) = dct;
            %    end
            %end
            
            Cr  = zeros(size(c));
            dCr = zeros(size(dc));
            for i = 1:3
                for m = 1:3
                    Rmis  = squeeze(R(m,i,:));                    
                    for n = 1:3
                        cmns = squeeze(c(m,n,:));
                        for j= 1:3
                            Rnjs  = squeeze(R(n,j,:));
                            Cij  = Rmis.*cmns.*Rnjs;
                            Cr(i,j,:)  = squeeze(Cr(i,j,:)) + Cij;
                            for ivar = 1:2
                                dCmns = squeeze(dc(m,n,ivar,:));
                                dCij = Rmis.*dCmns.*Rnjs;
                                dCr(i,j,ivar,:) = squeeze(dCr(i,j,ivar,:)) + dCij;
                            end
                        end
                    end
                end
            end
            obj.C = Cr;
            obj.dC = dCr;
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
        
        function computeDensity(obj,x)
            %rho = max(0.001,min(rho,0.999));
            %mx = max(0.01,min(sqrt(1-rho),0.99));
            %my = max(0.01,min(sqrt(1-rho),0.99));
            mx = x(:,1);
            my = x(:,2);
            [rho,drho] = obj.density.compute([mx,my]);
            obj.rho = rho;
            obj.drho = drho;
        end
        
    end
    
end