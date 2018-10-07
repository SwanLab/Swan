classdef SequentialLaminateHomogenizer < handle
    
    
    properties
        FiberConstitutiveTensor
        MatrixConstitutiveTensor
        ComplementaryTensor
    end
    
    methods
        
        function obj = SequentialLaminateHomogenizer()
            
            E_Fiber   = 1;
            E_Matrix  = 1e-3;
            nu_Fiber  = 1/3;
            nu_Matrix = 1/3;
            
            %E_Fiber = sym('E','positive');
            %nu_Fiber = sym('nu','positive');
            
            %mu = 2;%sym('mu','positive');
            %lambda = 1;%sym('lambda','positive');
            
            %E_Fiber = mu*(3*lambda+2*mu)/(lambda + mu);
            %nu_Fiber = lambda/(2*(lambda+mu));
            
            obj.FiberConstitutiveTensor  = IsotropicConstitutiveTensor3D(E_Fiber,nu_Fiber);
            obj.MatrixConstitutiveTensor = IsotropicConstitutiveTensor3D(E_Matrix,nu_Matrix);
            
            %d1 = sym('d1','real');
            %d2 = sym('d2','real');
            d1 = 1; d2 = 0; d3 = 1;
            direction = [d1 d2 d3];            
            direction = direction/norm(direction);
            
            obj.ComplementaryTensor = ComplementaryTensor(obj.FiberConstitutiveTensor,direction);
            obj.ComplementaryTensor.computeTensorVoigtInPlaneStress()
            
            %Check-isotropy
            C11 = obj.ComplementaryTensor.tensorVoigtInPlaneStress(1,1);
            C22 = obj.ComplementaryTensor.tensorVoigtInPlaneStress(2,2);
            C33 = obj.ComplementaryTensor.tensorVoigtInPlaneStress(3,3);
            C21 = obj.ComplementaryTensor.tensorVoigtInPlaneStress(2,1);
            simplify(C11-C22)
            simplify(C33 - 0.5*(C11+C22)-C21)
           
        end

        function sdf(obj)
            A =	m1 * ( (lambda+2*mu) - 1/mu*(lambda^2*e1y^2+(lambda+2*mu)^2*e1x^2) + K*((lambda+2*mu)*e1x^2+lambda*e1y^2)^2 ) ...
                +m2* ( (lambda+2*mu) - 1/mu*(lambda^2*e2y^2+(lambda+2*mu)^2*e2x^2) + K*((lambda+2*mu)*e2x^2+lambda*e2y^2)^2 );
            B =	m1 * ( (lambda+2*mu) - 1/mu*(lambda^2*e1x^2+(lambda+2*mu)^2*e1y^2) + K*((lambda+2*mu)*e1y^2+lambda*e1x^2)^2 ) ...
                +m2* ( (lambda+2*mu) - 1/mu*(lambda^2*e2x^2+(lambda+2*mu)^2*e2y^2) + K*((lambda+2*mu)*e2y^2+lambda*e2x^2)^2 );
            C =	m1 * ( 4*mu - 1/mu*(2*mu)^2 + K*(4*mu*e1x*e1y)^2 )/2 ...
                +m2* ( 4*mu - 1/mu*(2*mu)^2 + K*(4*mu*e2x*e2y)^2 )/2;
            D =  m1 * ( 2*lambda - 1/mu*(2*lambda*(lambda+2*mu)) + K*2*((lambda+2*mu)*e1y^2+lambda*e1x^2)*((lambda+2*mu)*e1x^2+lambda*e1y^2) )/2 ...
                +m2* ( 2*lambda - 1/mu*(2*lambda*(lambda+2*mu)) + K*2*((lambda+2*mu)*e2y^2+lambda*e2x^2)*((lambda+2*mu)*e2x^2+lambda*e2y^2) )/2;
            E =	m1 * ( -1/mu*(4*mu*(2*lambda+2*mu)*e1x*e1y) + K*2*(4*mu*e1x*e1y*((lambda+2*mu)*e1y^2+lambda*e1x^2)) )/(2*sqrt2) ...
                +m2* ( -1/mu*(4*mu*(2*lambda+2*mu)*e2x*e2y) + K*2*(4*mu*e2x*e2y*((lambda+2*mu)*e2y^2+lambda*e2x^2)) )/(2*sqrt2);
            F = m1 * ( -1/mu*(4*mu*(2*lambda+2*mu)*e1x*e1y) + K*2*(4*mu*e1x*e1y*((lambda+2*mu)*e1x^2+lambda*e1y^2)) )/(2*sqrt2) ...
                +m2* ( -1/mu*(4*mu*(2*lambda+2*mu)*e2x*e2y) + K*2*(4*mu*e2x*e2y*((lambda+2*mu)*e2x^2+lambda*e2y^2)) )/(2*sqrt2);
        end
        
    end

      
end

