classdef MixtureTheoryHomogenizer < handle
    
    properties (Access = public) 
        Ch
    end
    
    properties (Access = private)
        StiffTensor
        WeakTensor
        Theta
        Ex
        Ey
        nu_xy
        nu_yx
        mu        
    end
    
    methods (Access = public)
        
        function obj = MixtureTheoryHomogenizer(StiffTensor,WeakTensor,Theta)
            obj.init(StiffTensor,WeakTensor,Theta)
            obj.computeOrthotropicProperties()         
            obj.computeHomogenizedTensor()
        end
       
    end
    
    methods (Access = private)
        
        function init(obj,StiffTensor,WeakTensor,Theta)
            obj.StiffTensor = StiffTensor;
            obj.WeakTensor = WeakTensor;
            obj.Theta = Theta;            
        end
        
        function computeOrthotropicProperties(obj)
            E1 = obj.StiffTensor.E;
            E0 = obj.WeakTensor.E;
            nu1 = obj.StiffTensor.nu; 
            nu0 = obj.WeakTensor.nu; 
            mu1 = obj.StiffTensor.mu;
            mu0 = obj.WeakTensor.mu;
            
            obj.Ex = obj.serialize(E1,E0,obj.Theta);
            obj.Ey = obj.parelalize(E1,E0,obj.Theta);
            obj.nu_xy = obj.serialize(nu1,nu0,obj.Theta);           
            obj.nu_yx = obj.nu_xy*obj.Ey/obj.Ex;
            obj.mu = obj.parelalize(mu1,mu0,obj.Theta);  
        end
        
        function computeHomogenizedTensor(obj)
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
            obj.Ch  = C;
        end
        
        function fSerial = serialize(obj,f1,f0,rho)
            fSerial = rho*f1 + (1-rho)*f0;
        end
        
        function fParalel = parelalize(obj,f1,f0,rho)
            fParalel = 1/(rho/f1 + (1-rho)/f0);
        end
        
    end
    
    
end


