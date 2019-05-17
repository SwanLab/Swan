classdef HomogenizedVarComputerFromVademecum ...
        < HomogenizedVarComputer
    
    properties (Access = public)
        Ctensor
        density
    end
    
    properties (Access = private)
       Reps      
       Rot
       Cref
       dCref
    end
    
    methods (Access = public)
        
        function obj = HomogenizedVarComputerFromVademecum(cParams)
            s.fileName = [cParams.fileName,'Reduced'];
            v = VademecumVariablesLoader(s);
            v.load();
            obj.Ctensor = v.Ctensor;
            obj.density = v.density;
            obj.computeSymbolicRotationMatrix();
            obj.designVariable = cParams.designVariable;
        end
        
        function computeCtensor(obj,x)
            obj.obtainReferenceConsistutiveTensor(x);            
            obj.computeRotationMatrix();
            obj.rotateConstitutiveTensor();
        end
        
        function computeDensity(obj,x)
            mx = x(:,1);
            my = x(:,2);
            [rho,drho] = obj.density.compute([mx,my]);
            obj.rho = rho;
            obj.drho = drho;
        end        
        
    end
    
    methods (Access = private)
                
        function obtainReferenceConsistutiveTensor(obj,x)
            mx = x(:,1);
            my = x(:,2);
            [c,dc] = obj.Ctensor.compute([mx,my]);
            obj.Cref  = c;
            obj.dCref = permute(dc,[1 2 4 3]);            
        end        
        
        function computeRotationMatrix(obj)
            Rsym = obj.Reps;
            dir = obj.designVariable.alpha;            
            nx = squeeze(dir(1,:));
            ny = squeeze(dir(2,:));
            angle = squeeze(atan2(ny,nx));
            R = zeros(3,3,length(angle));
            for i = 1:3
                for j = 1:3
                    rij = matlabFunction(Rsym(i,j));
                    R(i,j,:) = rij(angle);
                end
            end
            obj.Rot = R;
        end        
        
        function rotateConstitutiveTensor(obj)
            cr  = obj.Cref;    
            dCr = obj.dCref;
            c  = zeros(size(cr));
            dc = zeros(size(dCr));
            nstre = size(cr,1);
            for i = 1:nstre
                for m = 1:nstre
                    Rmi  = obj.rotationComponent(m,i);  
                    for n = 1:nstre
                        Cmn  = squeeze(cr(m,n,:));
                        dCmn = squeeze(dCr(m,n,:,:));
                        for j = 1:nstre
                            Rnj = obj.rotationComponent(n,j);
                            Cij = Rmi.*Cmn.*Rnj;
                            c(i,j,:)  = squeeze(c(i,j,:)) + Cij;
                            for s = 1:2
                                dCmns(:,1) = dCmn(s,:);
                                dCij = Rmi.*dCmns.*Rnj;
                                dc(i,j,s,:) = squeeze(dc(i,j,s,:)) + dCij;
                            end
                        end
                    end
                end
            end      
            obj.C = c;
            obj.dC = dc;            
        end
        
        function Rij = rotationComponent(obj,i,j)
            R = obj.Rot;
            Rij = squeeze(R(i,j,:));            
        end        
        
        function a = computeSymbolicRotationMatrix(obj)
            S = StressPlaneStressVoigtTensor();
            alpha = sym('alpha','real');
            v = Vector3D();
            v.setValue([0 0 1]);
            factory = RotatorFactory();
            rotator = factory.create(S,alpha,v);
            a = rotator.rotationMatrix();  
            obj.Reps = simplify(inv(a)');
        end
        
    end
    
end