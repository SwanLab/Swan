classdef Rank2Homogenizer2D < handle
    
    properties
        Ch
    end
    
    properties (Access = private)
        E
        nu
        lambda2D
        mu
        lambda
        IsoTensor
        direction1
        direction2
        m1
        m2
        theta
    end
    
    methods
        
        function obj = Rank2Homogenizer2D(E,nu,dir1,dir2,m1,m2,theta,epsilon)
            obj.init(E,nu,dir1,dir2,m1,m2,theta);
            %Ca1 = obj.ComputeAnisotropicContribution(obj.direction1);
            Ca1 = obj.ComputeAnisotropicContribution2(obj.direction1);
            %Ca2 = obj.ComputeAnisotropicContribution(obj.direction2);
            Ca2 = obj.ComputeAnisotropicContribution2(obj.direction2);
            obj.computeCh(Ca1,Ca2,epsilon);
        end
        
    end
    
    methods (Access = private)
        function init(obj,E,nu,dir1,dir2,m1,m2,theta)
            obj.Ch = zeros(3,3);
            obj.E = E;
            obj.nu = nu;
            obj.direction1 = dir1;
            obj.direction2 = dir2; 
            obj.m1 = m1;
            obj.m2 = m2;
            obj.theta = theta;
            
            obj.IsoTensor = IsotropicConstitutiveTensor3D(E,nu);
            obj.computeLambda2DFromYoungAndPoisson();
            obj.computeMuFromYoungAndPoisson();
            obj.computeLambdaFromYoungAndPoisson();
        end
        
        function computeLambda2DFromYoungAndPoisson(obj)
            obj.lambda2D = obj.E*obj.nu/(1+obj.nu)/(1-obj.nu);
        end
        
        function computeMuFromYoungAndPoisson(obj)
            obj.mu = obj.E/(2*(1+obj.nu));           
        end
        
        function computeLambdaFromYoungAndPoisson(obj)
            obj.lambda = (obj.E*obj.nu/((1+obj.nu)*(1-2*obj.nu)));
        end

        
        function computeCh(obj,Cm1,Cm2,epsil)
           m1v = obj.m1;
           m2v = obj.m2;
           thet = obj.theta;
           
            E1 = obj.E;
            nu1 = obj.nu;
            E0 = epsil/(1-epsil)*obj.E;
            nu0 = obj.nu; 
           
           
            Tens1 = IsotropicConstitutiveTensor3D(E1,nu1);            
            Tens1.computeTensorVoigtInPlaneStress()
            Tens0  = IsotropicConstitutiveTensor3D(E0,nu0);            
            Tens0.computeTensorVoigtInPlaneStress()
           
            C1 = Tens1.tensorVoigtInPlaneStress();
            C0 = Tens0.tensorVoigtInPlaneStress();
            S0 = Inverter.invert(C0);
            S1 = Inverter.invert(C1);
           
            Cm = Cm1*m1v + Cm2*m2v;
            
            S01 = S0 - S1;
            C01 = Inverter.invert(S01);
            
            Ctheta = C01 + thet*Cm;   
            
            Stheta = Inverter.invert(Ctheta);
            
            Sh = S1 +(1-thet)*Stheta;            
            obj.Ch = Inverter.invert(Sh);
            
            
        end
        
        
        function CaVoigt = ComputeAnisotropicContribution(obj,dir)
            
            lam = obj.lambda;
            dim = 2;
            I = eye(dim);
            Em = dir(1:2)'*dir(1:2);
            
            c1 = 2*lam*obj.mu/(2*obj.mu+lam);
            c2 = 2*(lam+obj.mu)/lam;
            M1 = c1* [1 -1;
                      -1 c2];
            
            M2 = obj.mu*[1 -1
                  -1 0];
            
              Ca = zeros(dim,dim,dim,dim);
              
            for i = 1:dim
                for j = 1:dim
                    for k = 1:dim
                        for l = 1:dim
                            
                            d1 = I(i,j);
                            d2 = I(i,k);
                            d3 = I(i,l);
                            d4 = I(j,k);
                            d5 = I(j,l);
                            d6 = I(k,l);
                            
                            e1 = Em(i,j);
                            e2 = Em(i,k);
                            e3 = Em(i,l);
                            e4 = Em(j,k);
                            e5 = Em(j,l);
                            e6 = Em(k,l);
                            
                            
                            v1 = [d1 ;e1];
                            v2 = [d2 ;e2];
                            v3 = [d3 ;e3];
                            v4 = [d4 ;e4];
                            v5 = [d5 ;e5];
                            v6 = [d6 ;e6];
                            
                            
                            
                            t1 = v1'*M1*v6;
                            t2 = v2'*M2*v5;
                            t3 = v3'*M2*v4;
                            
                            Ca(i,j,k,l) = t1 + t2 + t3;
                        end
                    end
                end
            end
            

            
            CaVoigt = obj.computeVoigt(Ca);
            
            
        end
        
        function CaVoigt = computeVoigt(obj,Ca)
            CaVoigt = zeros(3,3);
            CaVoigt(1,1) = Ca(1,1,1,1);
            CaVoigt(1,2) = Ca(1,1,2,2);
            CaVoigt(2,2) = Ca(2,2,2,2);
            CaVoigt(1,3) = 0.5*Ca(1,1,1,2);
            CaVoigt(2,3) = 0.5*Ca(2,2,1,2);
            CaVoigt(3,3) = 0.5*Ca(1,2,1,2);
            
            
            CaVoigt(2,1) = CaVoigt(1,2);
            CaVoigt(3,1) = CaVoigt(1,3);
            CaVoigt(3,2) = CaVoigt(2,3);
        end
        
        function CaDir = ComputeAnisotropicContribution2(obj,dir)
            
            lam = obj.lambda;
            muT = obj.mu;
            dim = 2;
            K  = (lam+muT)/(2*muT+lam)/muT;
            I = eye(dim);
            dim = 2;
            Is = zeros(dim,dim,dim,dim);
            An1 = zeros(dim,dim,dim,dim);
            An2 = zeros(dim,dim,dim,dim);
            Ca = zeros(dim,dim,dim,dim);
            
            for i = 1:dim
                for j = 1:dim
                    for k = 1:dim
                        for l = 1:dim
                            Is(i,j,k,l) = lam*I(i,j)*I(k,l) + ...
                                +muT*(I(i,k)*I(j,l) + I(i,l)*I(j,k));
                        end
                    end
                end
            end
            
            
            
            for i = 1:dim
                for j = 1:dim
                    for k = 1:dim
                        for l = 1:dim

                            for p = 1:dim
                                for r = 1:dim
                                    for q = 1:dim
                                      An1(i,j,k,l) =  An1(i,j,k,l) - 1/muT*(Is(i,j,p,q)*dir(p)*Is(k,l,r,q)*dir(r));
                                    end
                                end
                            end
                            
                            for p = 1:dim
                                for q = 1:dim
                                    for r = 1:dim
                                        for s = 1:dim
                                            An2(i,j,k,l) =  An2(i,j,k,l) + K*(Is(p,q,k,l)*dir(p)*dir(q)*Is(r,s,i,j)*dir(r)*dir(s));
                                        end
                                    end
                                end
                            end
                            
                             Ca(i,j,k,l) = Is(i,j,k,l) + An1(i,j,k,l) + An2(i,j,k,l);
                        end
                    end
                end
            end
            
            
            CaDir = obj.computeVoigt(Ca); 
           
            
        end
        

        
        
    end
    
end
