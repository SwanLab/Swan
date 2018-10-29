classdef MixtureTheoryHomogenizer < handle
    
    properties (Access = public) 
        Ch
    end
    
    properties (Access = private)
        StiffTensor
        WeakTensor
        Angle
        Ex
        Ey
        nu_xy
        nu_yx
        mu        
        Direction
        ChHorizontal
        FractionVolume
    end
    
    methods (Access = public)
        
        function obj = MixtureTheoryHomogenizer(StiffTensor,WeakTensor,Dir,Angle,Vfrac)
            obj.init(StiffTensor,WeakTensor,Dir,Angle,Vfrac)
            obj.computeOrthotropicProperties() 
            obj.computeHorizontalHomogenizedTensor()
            obj.rotateHorizontalHomogenizedTensor()
        end
       
    end
    
    methods (Access = private)
        
        function init(obj,StiffTensor,WeakTensor,Dir,Angle,Vfrac)
            obj.StiffTensor = StiffTensor;
            obj.WeakTensor = WeakTensor;
            obj.Angle = Angle;            
            obj.Direction = Dir;
            obj.FractionVolume = Vfrac;
        end
        
        function computeOrthotropicProperties(obj)
            E1 = obj.StiffTensor.E;
            E0 = obj.WeakTensor.E;
            nu1 = obj.StiffTensor.nu; 
            nu0 = obj.WeakTensor.nu; 
            mu1 = obj.StiffTensor.mu;
            mu0 = obj.WeakTensor.mu;
            Vfrac = obj.FractionVolume;
            
            
            obj.Ex = obj.serialize(E1,E0,Vfrac);
            obj.Ey = obj.parelalize(E1,E0,Vfrac);
            obj.nu_xy = obj.serialize(nu1,nu0,Vfrac);           
            obj.nu_yx = obj.nu_xy*obj.Ey/obj.Ex;
            obj.mu = obj.parelalize(mu1,mu0,Vfrac);
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
        
        function fSerial = serialize(obj,f1,f0,rho)
            fSerial = rho*f1 + (1-rho)*f0;
        end
        
        function fParalel = parelalize(obj,f1,f0,rho)
            fParalel = 1/(rho/f1 + (1-rho)/f0);
        end
        
        function rotateHorizontalHomogenizedTensor(obj)
            d = obj.Direction;
            a = obj.Angle;
            ChHor = obj.ChHorizontal;
%             rotator = StressVoigtPlaneStressRotator(a,d);
%             R = rotator.getRotationMatrix();            
%             obj.Ch = R*ChHor*(R');     
            
            C = FourthOrderVoigtTensor();
            C.setValue(ChHor);            
            obj.Ch = Rotator.rotate(C,a,d);
        end
        
    end
    
    
end


