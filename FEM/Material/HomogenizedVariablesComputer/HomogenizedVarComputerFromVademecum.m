classdef HomogenizedVarComputerFromVademecum ...
        < HomogenizedVarComputer
    
    properties (Access = public)
        Ctensor
        density
    end
    
    properties (Access = private)
       Rsymbolic
       RsymbolicEps       
    end
    
    methods (Access = public)
        
        function obj = HomogenizedVarComputerFromVademecum(cParams)
            s.fileName = [cParams.fileName,'Reduced'];
            v = VademecumVariablesLoader(s);
            v.load();
            obj.Ctensor = v.Ctensor;
            obj.density = v.density;
            obj.createSymbolicRotationMatrix();
            obj.designVariable = cParams.designVariable;
        end
        
        function R = computeRotationMatrix2(obj)
           b = obj.designVariable.beta;            
           b1 = -b(1,:);
           b2 = b(2,:);
           R(1,1,:) = (1 - b1)/2;
           R(1,2,:) = (1 + b1)/2;
           R(1,3,:) = -b2;
           R(2,1,:) = (1 + b1)/2;
           R(2,2,:) = (1 - b1)/2;
           R(2,3,:) = b2;
           R(3,1,:) = b2/2;
           R(3,2,:) = -b2/2;
           R(3,3,:) = -b1;           
        end
        
        function computeCtensor(obj,x)
            
            Rs = obj.computeRotatorMatrix(obj.Rsymbolic);
            Rp = obj.computeRotationMatrix2();
            
            Re = obj.computeRotatorMatrix(obj.RsymbolicEps);
            
            R = Re;
           % Rs = Rs2;

            mx = x(:,1);
            my = x(:,2);
            [c,dc] = obj.Ctensor.compute([mx,my]);
            dc = permute(dc,[1 2 4 3]);

            
            Cr  = zeros(size(c));
            dCr = zeros(size(dc));
            for i = 1:3
                for m = 1:3
                    Rmi  = squeeze(R(m,i,:));  
                    Rim = squeeze(R(i,m,:));                                        
                    for n = 1:3
                        Cmn = squeeze(c(m,n,:));
                        for j= 1:3
                            Rnj  = squeeze(R(n,j,:));
                            Rjn  = squeeze(R(j,n,:));
                            
                            Cij  = Rmi.*Cmn.*Rnj;
                            %Cij  = Rim.*Cmn.*Rjn;
                            
                            Cr(i,j,:)  = squeeze(Cr(i,j,:)) + Cij;
                            for ivar = 1:2
                                dCmn = squeeze(dc(m,n,ivar,:));
                                dCij = Rmi.*dCmn.*Rnj;
                                %dCij = Rim.*dCmn.*Rjn;
                                dCr(i,j,ivar,:) = squeeze(dCr(i,j,ivar,:)) + dCij;
                            end
                        end
                    end
                end
            end
            
            obj.C = Cr;
            obj.dC = dCr;
        end
        
        function R = computeRotatorMatrix(obj,Rs)
            dir = obj.designVariable.alpha;            
            nx = squeeze(dir(1,:));
            ny = squeeze(dir(2,:));
            angle = squeeze(atan2(ny,nx));
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
            alpha = sym('alpha','real');
            v = Vector3D();
            v.setValue([0 0 1]);
            factory = RotatorFactory();
            rotator = factory.create(S,alpha,v);
            a = rotator.rotationMatrix();  
            obj.Rsymbolic = a;
            obj.RsymbolicEps = simplify(inv(a)');
        end
        
        function computeDensity(obj,x)
            mx = x(:,1);
            my = x(:,2);
            [rho,drho] = obj.density.compute([mx,my]);
            obj.rho = rho;
            obj.drho = drho;
        end
        
    end
    
end