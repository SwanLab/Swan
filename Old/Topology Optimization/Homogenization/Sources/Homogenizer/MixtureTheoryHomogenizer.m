classdef MixtureTheoryHomogenizer < handle
    
    properties (Access = public) 
        Ch
    end
    
    properties (Access = private)
        C1
        C0
        angle
        Ex
        Ey
        nu_xy
        nu_yx
        mu
        dir
        ChHorizontal
        Vfrac
    end
    
    methods (Access = public)
        
        function obj = MixtureTheoryHomogenizer(C1,C0,Dir,angle,Vfrac)
            obj.init(C1,C0,Dir,angle,Vfrac)
            obj.computeOrthotropicProperties() 
            obj.computeHorizontalHomogenizedTensor()
            obj.rotateHorizontalHomogenizedTensor()
        end
       
    end
    
    methods (Access = private)
        
        function init(obj,C1,C0,Dir,angle,Vfrac)
            obj.C1 = C1;
            obj.C0 = C0;
            obj.angle = angle;
            obj.dir = Dir;
            obj.Vfrac = Vfrac;
        end
        
        function computeOrthotropicProperties(obj)
            E1 = obj.C1.getYoung();
            E0 = obj.C0.getYoung();
            nu1 = obj.C1.getPoisson();
            nu0 = obj.C0.getPoisson();
            mu1 = obj.C1.getMu();
            mu0 = obj.C0.getMu();
            Vol = obj.Vfrac;
            
            
            obj.Ex = obj.serialize(E1,E0,Vol);
            obj.Ey = obj.parelalize(E1,E0,Vol);
            obj.nu_xy = obj.serialize(nu1,nu0,Vol);
            obj.nu_yx = obj.nu_xy*obj.Ey/obj.Ex;
            obj.mu = obj.parelalize(mu1,mu0,Vol);
        end

        function computeHorizontalHomogenizedTensor(obj)
            E1    = obj.Ex;
            E2    = obj.Ey;
            nu_12 = obj.nu_xy;
            nu_21 = obj.nu_yx;
            Mu    = obj.mu;
            
            C = zeros(3,3);
            C(1,1) = E1/(1-nu_12*nu_21);
            C(2,2) = E2/(1-nu_12*nu_21);
            C(1,2) = E1*nu_21/(1-nu_12*nu_21);
            C(2,1) = E2*nu_12/(1-nu_12*nu_21);
            C(3,3) = Mu;
            obj.ChHorizontal  = C;
        end
        
        function rotateHorizontalHomogenizedTensor(obj)
            d = obj.dir;
            a = obj.angle;
            ChHor = obj.ChHorizontal;
            C = SymmetricFourthOrderPlaneStressVoigtTensor();
            C.setValue(ChHor);
            obj.Ch = Rotator.rotate(C,a,d);
        end

    end
        
    methods (Access = private, Static)
        function fSerial = serialize(f1,f0,rho)
            fSerial = rho*f1 + (1-rho)*f0;
        end
        
        function fParalel = parelalize(f1,f0,rho)
            fParalel = 1/(rho/f1 + (1-rho)/f0);
        end
        
    end
    
end
