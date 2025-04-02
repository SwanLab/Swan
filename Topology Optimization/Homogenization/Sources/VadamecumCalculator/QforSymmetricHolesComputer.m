classdef QforSymmetricHolesComputer < handle
    
    properties (Access = private)
        A
        B        
        C
        D
        E
    end
    
    properties (Access = private)
        rho0
        rho1
        rhoM
        rhoS
        qRho0
        qRho1
        qRhoM
        qRhoS
    end
    
    properties (Access = private)
        
    end
    
    methods (Access = public)
        
        function obj = QforSymmetricHolesComputer(cParams)
            obj.init(cParams)            
        end
        
        function qMin = compute(obj)
            obj.computeCoefficients();
            qMin = @(rho) (obj.A*rho.^2 + obj.B*rho + obj.C)./(obj.D*rho.^2 + obj.E*rho + 1);
        end        
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.rho0  = cParams.rho0;
            obj.rho1  = cParams.rho1;
            obj.rhoM  = cParams.rhoM;
            obj.rhoS  = cParams.rhoS;
            obj.qRho0 = cParams.qRho0;
            obj.qRho1 = cParams.qRho1;
            obj.qRhoM = cParams.qRhoM;
            obj.qRhoS = cParams.qRhoS;
        end
        
        function computeCoefficients(obj)
            r0 = obj.rho0;
            r1 = obj.rho1;
            rS = obj.rhoM;
            rM = obj.rhoS;
            q0 = obj.qRho0;
            q1 = obj.qRho1;
            qM = obj.qRhoM;
            qS = obj.qRhoS;            
            obj.A = -(q0*qM^2*r0^2*r1 - q1*qM^2*r0*r1^2 + q0*qM^2*r1*rS^2 - q1*qM^2*r0*rS^2 - q0*qM^2*r0^2*rM + q1*qM^2*r1^2*rM - q0*qM^2*rS^2*rM + q1*qM^2*rS^2*rM + qM^2*qS*r0*rS^2 - qM^2*qS*r1*rS^2 + qM^2*qS*r0*rM^2 - qM^2*qS*r1*rM^2 + q0*q1*qM*r0*r1^2 - q0*q1*qM*r0^2*r1 + q0*q1*qM*r0*rS^2 - q0*q1*qM*r1*rS^2 + q0*q1*qM*r0^2*rM - q0*q1*qM*r1^2*rM - q0*q1*qS*r0*r1^2 + q0*q1*qS*r0^2*r1 + q0*q1*qS*r0*rM^2 - q0*q1*qS*r0^2*rM - q0*q1*qS*r1*rM^2 + q0*q1*qS*r1^2*rM - q0*qM*qS*r0^2*r1 + q1*qM*qS*r0*r1^2 - q0*qM*qS*r0*rS^2 + q1*qM*qS*r1*rS^2 - q0*qM*qS*r0*rM^2 + q0*qM*qS*r0^2*rM + q0*qM*qS*r1*rM^2 - q1*qM*qS*r0*rM^2 + q1*qM*qS*r1*rM^2 - q1*qM*qS*r1^2*rM + q0*qM*qS*rS^2*rM - q1*qM*qS*rS^2*rM - 2*q0*qM^2*r0*r1*rS + 2*q1*qM^2*r0*r1*rS + 2*q0*qM^2*r0*rS*rM - 2*q1*qM^2*r1*rS*rM - 2*qM^2*qS*r0*rS*rM + 2*qM^2*qS*r1*rS*rM - 2*q0*q1*qM*r0*rS*rM + 2*q0*q1*qM*r1*rS*rM + 2*q0*qM*qS*r0*r1*rS - 2*q1*qM*qS*r0*r1*rS - 2*q0*qM*qS*r1*rS*rM + 2*q1*qM*qS*r0*rS*rM)/(qM^2*r0*r1^2*rS^2 - qM^2*r0^2*r1*rS^2 - qM^2*r0*rS^2*rM^2 + qM^2*r0^2*rS^2*rM + qM^2*r1*rS^2*rM^2 - qM^2*r1^2*rS^2*rM + q0*q1*r0*r1^2*rS^2 - q0*q1*r0^2*r1*rS^2 + q0*q1*r0*r1^2*rM^2 - q0*q1*r0^2*r1*rM^2 - q0*qM*r0*r1^2*rS^2 + q0*qM*r0^2*r1*rS^2 - q1*qM*r0*r1^2*rS^2 + q1*qM*r0^2*r1*rS^2 + q0*qM*r0^2*r1*rM^2 - q0*qM*r0^2*r1^2*rM - q1*qM*r0*r1^2*rM^2 + q1*qM*r0^2*r1^2*rM + q0*qM*r0*rS^2*rM^2 - q0*qM*r0^2*rS^2*rM - q1*qM*r1*rS^2*rM^2 + q1*qM*r1^2*rS^2*rM - q0*qS*r0*r1^2*rM^2 + q0*qS*r0^2*r1^2*rM + q1*qS*r0^2*r1*rM^2 - q1*qS*r0^2*r1^2*rM - q0*qS*r0*rS^2*rM^2 + q0*qS*r0^2*rS^2*rM + q1*qS*r1*rS^2*rM^2 - q1*qS*r1^2*rS^2*rM + qM*qS*r0*r1^2*rM^2 - qM*qS*r0^2*r1*rM^2 + qM*qS*r0*rS^2*rM^2 - qM*qS*r0^2*rS^2*rM - qM*qS*r1*rS^2*rM^2 + qM*qS*r1^2*rS^2*rM - 2*q0*q1*r0*r1^2*rS*rM + 2*q0*q1*r0^2*r1*rS*rM - 2*q0*qM*r0*r1*rS*rM^2 + 2*q0*qM*r0*r1^2*rS*rM + 2*q1*qM*r0*r1*rS*rM^2 - 2*q1*qM*r0^2*r1*rS*rM + 2*q0*qS*r0*r1*rS*rM^2 - 2*q0*qS*r0^2*r1*rS*rM - 2*q1*qS*r0*r1*rS*rM^2 + 2*q1*qS*r0*r1^2*rS*rM - 2*qM*qS*r0*r1^2*rS*rM + 2*qM*qS*r0^2*r1*rS*rM);
            obj.B = (q0*qM^2*r0^2*r1^2 - q1*qM^2*r0^2*r1^2 + q0*qM^2*r1^2*rS^2 - q1*qM^2*r0^2*rS^2 - q0*qM^2*r0^2*rM^2 + q1*qM^2*r1^2*rM^2 - q0*qM^2*rS^2*rM^2 + q1*qM^2*rS^2*rM^2 + qM^2*qS*r0^2*rS^2 - qM^2*qS*r1^2*rS^2 + qM^2*qS*r0^2*rM^2 - qM^2*qS*r1^2*rM^2 + q0*q1*qM*r0^2*rS^2 - q0*q1*qM*r1^2*rS^2 + q0*q1*qM*r0^2*rM^2 - q0*q1*qM*r1^2*rM^2 - q0*qM*qS*r0^2*r1^2 + q1*qM*qS*r0^2*r1^2 - q0*qM*qS*r0^2*rS^2 + q1*qM*qS*r1^2*rS^2 + q0*qM*qS*r1^2*rM^2 - q1*qM*qS*r0^2*rM^2 + q0*qM*qS*rS^2*rM^2 - q1*qM*qS*rS^2*rM^2 - 2*q0*qM^2*r0*r1^2*rS + 2*q1*qM^2*r0^2*r1*rS + 2*q0*qM^2*r0*rS*rM^2 - 2*q1*qM^2*r1*rS*rM^2 - 2*qM^2*qS*r0^2*rS*rM + 2*qM^2*qS*r1^2*rS*rM + 2*q0*q1*qM*r0*r1^2*rS - 2*q0*q1*qM*r0^2*r1*rS - 2*q0*q1*qM*r0*rS*rM^2 + 2*q0*q1*qM*r1*rS*rM^2 - 2*q0*q1*qS*r0*r1^2*rS + 2*q0*q1*qS*r0^2*r1*rS + 2*q0*q1*qS*r0*rS*rM^2 - 2*q0*q1*qS*r0^2*rS*rM - 2*q0*q1*qS*r1*rS*rM^2 + 2*q0*q1*qS*r1^2*rS*rM + 2*q0*qM*qS*r0*r1^2*rS - 2*q1*qM*qS*r0^2*r1*rS - 2*q0*qM*qS*r0*rS*rM^2 + 2*q0*qM*qS*r0^2*rS*rM - 2*q0*qM*qS*r1^2*rS*rM + 2*q1*qM*qS*r0^2*rS*rM + 2*q1*qM*qS*r1*rS*rM^2 - 2*q1*qM*qS*r1^2*rS*rM)/(qM^2*r0*r1^2*rS^2 - qM^2*r0^2*r1*rS^2 - qM^2*r0*rS^2*rM^2 + qM^2*r0^2*rS^2*rM + qM^2*r1*rS^2*rM^2 - qM^2*r1^2*rS^2*rM + q0*q1*r0*r1^2*rS^2 - q0*q1*r0^2*r1*rS^2 + q0*q1*r0*r1^2*rM^2 - q0*q1*r0^2*r1*rM^2 - q0*qM*r0*r1^2*rS^2 + q0*qM*r0^2*r1*rS^2 - q1*qM*r0*r1^2*rS^2 + q1*qM*r0^2*r1*rS^2 + q0*qM*r0^2*r1*rM^2 - q0*qM*r0^2*r1^2*rM - q1*qM*r0*r1^2*rM^2 + q1*qM*r0^2*r1^2*rM + q0*qM*r0*rS^2*rM^2 - q0*qM*r0^2*rS^2*rM - q1*qM*r1*rS^2*rM^2 + q1*qM*r1^2*rS^2*rM - q0*qS*r0*r1^2*rM^2 + q0*qS*r0^2*r1^2*rM + q1*qS*r0^2*r1*rM^2 - q1*qS*r0^2*r1^2*rM - q0*qS*r0*rS^2*rM^2 + q0*qS*r0^2*rS^2*rM + q1*qS*r1*rS^2*rM^2 - q1*qS*r1^2*rS^2*rM + qM*qS*r0*r1^2*rM^2 - qM*qS*r0^2*r1*rM^2 + qM*qS*r0*rS^2*rM^2 - qM*qS*r0^2*rS^2*rM - qM*qS*r1*rS^2*rM^2 + qM*qS*r1^2*rS^2*rM - 2*q0*q1*r0*r1^2*rS*rM + 2*q0*q1*r0^2*r1*rS*rM - 2*q0*qM*r0*r1*rS*rM^2 + 2*q0*qM*r0*r1^2*rS*rM + 2*q1*qM*r0*r1*rS*rM^2 - 2*q1*qM*r0^2*r1*rS*rM + 2*q0*qS*r0*r1*rS*rM^2 - 2*q0*qS*r0^2*r1*rS*rM - 2*q1*qS*r0*r1*rS*rM^2 + 2*q1*qS*r0*r1^2*rS*rM - 2*qM*qS*r0*r1^2*rS*rM + 2*qM*qS*r0^2*r1*rS*rM);
            obj.C = (q0*qM^2*r0^2*r1*rM^2 - q0*qM^2*r0^2*r1^2*rM - q1*qM^2*r0*r1^2*rM^2 + q1*qM^2*r0^2*r1^2*rM + q0*qM^2*r1*rS^2*rM^2 - q0*qM^2*r1^2*rS^2*rM - q1*qM^2*r0*rS^2*rM^2 + q1*qM^2*r0^2*rS^2*rM + qM^2*qS*r0*r1^2*rS^2 - qM^2*qS*r0^2*r1*rS^2 + qM^2*qS*r0*r1^2*rM^2 - qM^2*qS*r0^2*r1*rM^2 + q0*q1*qM*r0*r1^2*rM^2 - q0*q1*qM*r0^2*r1*rM^2 + q0*q1*qM*r0*rS^2*rM^2 - q0*q1*qM*r0^2*rS^2*rM - q0*q1*qM*r1*rS^2*rM^2 + q0*q1*qM*r1^2*rS^2*rM + q0*q1*qS*r0*r1^2*rS^2 - q0*q1*qS*r0^2*r1*rS^2 - q0*q1*qS*r0*rS^2*rM^2 + q0*q1*qS*r0^2*rS^2*rM + q0*q1*qS*r1*rS^2*rM^2 - q0*q1*qS*r1^2*rS^2*rM - q0*qM*qS*r0*r1^2*rS^2 + q0*qM*qS*r0^2*r1*rS^2 - q1*qM*qS*r0*r1^2*rS^2 + q1*qM*qS*r0^2*r1*rS^2 - q0*qM*qS*r0*r1^2*rM^2 + q0*qM*qS*r0^2*r1^2*rM + q1*qM*qS*r0^2*r1*rM^2 - q1*qM*qS*r0^2*r1^2*rM - q0*qM*qS*r1*rS^2*rM^2 + q0*qM*qS*r1^2*rS^2*rM + q1*qM*qS*r0*rS^2*rM^2 - q1*qM*qS*r0^2*rS^2*rM - 2*q0*qM^2*r0*r1*rS*rM^2 + 2*q0*qM^2*r0*r1^2*rS*rM + 2*q1*qM^2*r0*r1*rS*rM^2 - 2*q1*qM^2*r0^2*r1*rS*rM - 2*qM^2*qS*r0*r1^2*rS*rM + 2*qM^2*qS*r0^2*r1*rS*rM - 2*q0*q1*qM*r0*r1^2*rS*rM + 2*q0*q1*qM*r0^2*r1*rS*rM + 2*q0*qM*qS*r0*r1*rS*rM^2 - 2*q0*qM*qS*r0^2*r1*rS*rM - 2*q1*qM*qS*r0*r1*rS*rM^2 + 2*q1*qM*qS*r0*r1^2*rS*rM)/(qM^2*r0*r1^2*rS^2 - qM^2*r0^2*r1*rS^2 - qM^2*r0*rS^2*rM^2 + qM^2*r0^2*rS^2*rM + qM^2*r1*rS^2*rM^2 - qM^2*r1^2*rS^2*rM + q0*q1*r0*r1^2*rS^2 - q0*q1*r0^2*r1*rS^2 + q0*q1*r0*r1^2*rM^2 - q0*q1*r0^2*r1*rM^2 - q0*qM*r0*r1^2*rS^2 + q0*qM*r0^2*r1*rS^2 - q1*qM*r0*r1^2*rS^2 + q1*qM*r0^2*r1*rS^2 + q0*qM*r0^2*r1*rM^2 - q0*qM*r0^2*r1^2*rM - q1*qM*r0*r1^2*rM^2 + q1*qM*r0^2*r1^2*rM + q0*qM*r0*rS^2*rM^2 - q0*qM*r0^2*rS^2*rM - q1*qM*r1*rS^2*rM^2 + q1*qM*r1^2*rS^2*rM - q0*qS*r0*r1^2*rM^2 + q0*qS*r0^2*r1^2*rM + q1*qS*r0^2*r1*rM^2 - q1*qS*r0^2*r1^2*rM - q0*qS*r0*rS^2*rM^2 + q0*qS*r0^2*rS^2*rM + q1*qS*r1*rS^2*rM^2 - q1*qS*r1^2*rS^2*rM + qM*qS*r0*r1^2*rM^2 - qM*qS*r0^2*r1*rM^2 + qM*qS*r0*rS^2*rM^2 - qM*qS*r0^2*rS^2*rM - qM*qS*r1*rS^2*rM^2 + qM*qS*r1^2*rS^2*rM - 2*q0*q1*r0*r1^2*rS*rM + 2*q0*q1*r0^2*r1*rS*rM - 2*q0*qM*r0*r1*rS*rM^2 + 2*q0*qM*r0*r1^2*rS*rM + 2*q1*qM*r0*r1*rS*rM^2 - 2*q1*qM*r0^2*r1*rS*rM + 2*q0*qS*r0*r1*rS*rM^2 - 2*q0*qS*r0^2*r1*rS*rM - 2*q1*qS*r0*r1*rS*rM^2 + 2*q1*qS*r0*r1^2*rS*rM - 2*qM*qS*r0*r1^2*rS*rM + 2*qM*qS*r0^2*r1*rS*rM);
            obj.D = (qM^2*r0*r1^2 - qM^2*r0^2*r1 - qM^2*r0*rM^2 + qM^2*r0^2*rM + qM^2*r1*rM^2 - qM^2*r1^2*rM - q0*q1*r0*rS^2 + q0*q1*r1*rS^2 - q0*q1*r0*rM^2 + q0*q1*r1*rM^2 - q0*qM*r0*r1^2 + q1*qM*r0^2*r1 - q0*qM*r1*rS^2 + q1*qM*r0*rS^2 + q0*qM*r0*rM^2 - q0*qM*r1*rM^2 + q0*qM*r1^2*rM + q1*qM*r0*rM^2 - q1*qM*r0^2*rM - q1*qM*r1*rM^2 + q0*qM*rS^2*rM - q1*qM*rS^2*rM + q0*qS*r0*r1^2 - q1*qS*r0^2*r1 + q0*qS*r0*rS^2 - q1*qS*r1*rS^2 - q0*qS*r1^2*rM + q1*qS*r0^2*rM - q0*qS*rS^2*rM + q1*qS*rS^2*rM - qM*qS*r0*r1^2 + qM*qS*r0^2*r1 - qM*qS*r0*rS^2 + qM*qS*r1*rS^2 - qM*qS*r0^2*rM + qM*qS*r1^2*rM + 2*q0*q1*r0*rS*rM - 2*q0*q1*r1*rS*rM + 2*q0*qM*r0*r1*rS - 2*q1*qM*r0*r1*rS - 2*q0*qM*r0*rS*rM + 2*q1*qM*r1*rS*rM - 2*q0*qS*r0*r1*rS + 2*q1*qS*r0*r1*rS + 2*q0*qS*r1*rS*rM - 2*q1*qS*r0*rS*rM + 2*qM*qS*r0*rS*rM - 2*qM*qS*r1*rS*rM)/(qM^2*r0*r1^2*rS^2 - qM^2*r0^2*r1*rS^2 - qM^2*r0*rS^2*rM^2 + qM^2*r0^2*rS^2*rM + qM^2*r1*rS^2*rM^2 - qM^2*r1^2*rS^2*rM + q0*q1*r0*r1^2*rS^2 - q0*q1*r0^2*r1*rS^2 + q0*q1*r0*r1^2*rM^2 - q0*q1*r0^2*r1*rM^2 - q0*qM*r0*r1^2*rS^2 + q0*qM*r0^2*r1*rS^2 - q1*qM*r0*r1^2*rS^2 + q1*qM*r0^2*r1*rS^2 + q0*qM*r0^2*r1*rM^2 - q0*qM*r0^2*r1^2*rM - q1*qM*r0*r1^2*rM^2 + q1*qM*r0^2*r1^2*rM + q0*qM*r0*rS^2*rM^2 - q0*qM*r0^2*rS^2*rM - q1*qM*r1*rS^2*rM^2 + q1*qM*r1^2*rS^2*rM - q0*qS*r0*r1^2*rM^2 + q0*qS*r0^2*r1^2*rM + q1*qS*r0^2*r1*rM^2 - q1*qS*r0^2*r1^2*rM - q0*qS*r0*rS^2*rM^2 + q0*qS*r0^2*rS^2*rM + q1*qS*r1*rS^2*rM^2 - q1*qS*r1^2*rS^2*rM + qM*qS*r0*r1^2*rM^2 - qM*qS*r0^2*r1*rM^2 + qM*qS*r0*rS^2*rM^2 - qM*qS*r0^2*rS^2*rM - qM*qS*r1*rS^2*rM^2 + qM*qS*r1^2*rS^2*rM - 2*q0*q1*r0*r1^2*rS*rM + 2*q0*q1*r0^2*r1*rS*rM - 2*q0*qM*r0*r1*rS*rM^2 + 2*q0*qM*r0*r1^2*rS*rM + 2*q1*qM*r0*r1*rS*rM^2 - 2*q1*qM*r0^2*r1*rS*rM + 2*q0*qS*r0*r1*rS*rM^2 - 2*q0*qS*r0^2*r1*rS*rM - 2*q1*qS*r0*r1*rS*rM^2 + 2*q1*qS*r0*r1^2*rS*rM - 2*qM*qS*r0*r1^2*rS*rM + 2*qM*qS*r0^2*r1*rS*rM);
            obj.E = (q0*q1*r0^2*rS^2 - q0*q1*r1^2*rS^2 + q0*q1*r0^2*rM^2 - q0*q1*r1^2*rM^2 + q0*qM*r0^2*r1^2 - q1*qM*r0^2*r1^2 + q0*qM*r1^2*rS^2 - q1*qM*r0^2*rS^2 - q0*qM*r0^2*rM^2 + q1*qM*r1^2*rM^2 - q0*qM*rS^2*rM^2 + q1*qM*rS^2*rM^2 - q0*qS*r0^2*r1^2 + q1*qS*r0^2*r1^2 - q0*qS*r0^2*rS^2 + q1*qS*r1^2*rS^2 + q0*qS*r1^2*rM^2 - q1*qS*r0^2*rM^2 + q0*qS*rS^2*rM^2 - q1*qS*rS^2*rM^2 + qM*qS*r0^2*rS^2 - qM*qS*r1^2*rS^2 + qM*qS*r0^2*rM^2 - qM*qS*r1^2*rM^2 - 2*qM^2*r0*r1^2*rS + 2*qM^2*r0^2*r1*rS + 2*qM^2*r0*rS*rM^2 - 2*qM^2*r0^2*rS*rM - 2*qM^2*r1*rS*rM^2 + 2*qM^2*r1^2*rS*rM - 2*q0*q1*r0^2*rS*rM + 2*q0*q1*r1^2*rS*rM - 2*q0*qM*r0^2*r1*rS + 2*q1*qM*r0*r1^2*rS + 2*q0*qM*r0^2*rS*rM + 2*q0*qM*r1*rS*rM^2 - 2*q0*qM*r1^2*rS*rM - 2*q1*qM*r0*rS*rM^2 + 2*q1*qM*r0^2*rS*rM - 2*q1*qM*r1^2*rS*rM + 2*q0*qS*r0^2*r1*rS - 2*q1*qS*r0*r1^2*rS - 2*q0*qS*r1*rS*rM^2 + 2*q1*qS*r0*rS*rM^2 + 2*qM*qS*r0*r1^2*rS - 2*qM*qS*r0^2*r1*rS - 2*qM*qS*r0*rS*rM^2 + 2*qM*qS*r1*rS*rM^2)/(qM^2*r0*r1^2*rS^2 - qM^2*r0^2*r1*rS^2 - qM^2*r0*rS^2*rM^2 + qM^2*r0^2*rS^2*rM + qM^2*r1*rS^2*rM^2 - qM^2*r1^2*rS^2*rM + q0*q1*r0*r1^2*rS^2 - q0*q1*r0^2*r1*rS^2 + q0*q1*r0*r1^2*rM^2 - q0*q1*r0^2*r1*rM^2 - q0*qM*r0*r1^2*rS^2 + q0*qM*r0^2*r1*rS^2 - q1*qM*r0*r1^2*rS^2 + q1*qM*r0^2*r1*rS^2 + q0*qM*r0^2*r1*rM^2 - q0*qM*r0^2*r1^2*rM - q1*qM*r0*r1^2*rM^2 + q1*qM*r0^2*r1^2*rM + q0*qM*r0*rS^2*rM^2 - q0*qM*r0^2*rS^2*rM - q1*qM*r1*rS^2*rM^2 + q1*qM*r1^2*rS^2*rM - q0*qS*r0*r1^2*rM^2 + q0*qS*r0^2*r1^2*rM + q1*qS*r0^2*r1*rM^2 - q1*qS*r0^2*r1^2*rM - q0*qS*r0*rS^2*rM^2 + q0*qS*r0^2*rS^2*rM + q1*qS*r1*rS^2*rM^2 - q1*qS*r1^2*rS^2*rM + qM*qS*r0*r1^2*rM^2 - qM*qS*r0^2*r1*rM^2 + qM*qS*r0*rS^2*rM^2 - qM*qS*r0^2*rS^2*rM - qM*qS*r1*rS^2*rM^2 + qM*qS*r1^2*rS^2*rM - 2*q0*q1*r0*r1^2*rS*rM + 2*q0*q1*r0^2*r1*rS*rM - 2*q0*qM*r0*r1*rS*rM^2 + 2*q0*qM*r0*r1^2*rS*rM + 2*q1*qM*r0*r1*rS*rM^2 - 2*q1*qM*r0^2*r1*rS*rM + 2*q0*qS*r0*r1*rS*rM^2 - 2*q0*qS*r0^2*r1*rS*rM - 2*q1*qS*r0*r1*rS*rM^2 + 2*q1*qS*r0*r1^2*rS*rM - 2*qM*qS*r0*r1^2*rS*rM + 2*qM*qS*r0^2*r1*rS*rM);
        end
        
    end
    
end