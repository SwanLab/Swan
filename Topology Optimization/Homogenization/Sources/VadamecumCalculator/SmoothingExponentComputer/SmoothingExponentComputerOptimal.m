classdef SmoothingExponentComputerOptimal < SmoothingExponentComputer
    
    properties (Access = private)
        rho0
        rhoS        
        rhoM
        rho1        
        qRho0
        qRhoS        
        qRhoM          
        qRho1
        m1Max
        m2Max
        alpha
        beta
    end
    
    properties (Access = private)
        m1
        m2        
    end
    
    methods (Access = public)
        
        function obj = SmoothingExponentComputerOptimal(cParams)
            obj.init(cParams);
        end
        
        function setParamsValues(obj,alpha,beta,rM,qS)
            obj.beta  = beta;
            obj.alpha = alpha;
            obj.rhoM  = rM;
            obj.qRhoS = qS;
        end
        
        function q = computeQ(obj,xi,rho)
            qSym = obj.computeQForSymmetricHoles(rho);
            qElo = obj.computeQForElongatedHoles(rho);
            xiS  = obj.computeXiScaled(xi,rho);
            q = qSym + xiS.^(obj.alpha).*(qElo - qSym);
            q = max(min(obj.qRho0,q),obj.qRhoM);
        end
        
    end
    
    methods (Access = protected)
        
        function computeExponent(obj)
            obj.computeExponentAnalytical();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.m1 = cParams.m1;
            obj.m2 = cParams.m2;
            obj.computeXi();            
            obj.loadParameters();
        end
        
        function loadParameters(obj)
            obj.m1Max  = 0.99;
            obj.m2Max  = 0.99;
            obj.rho0   = 0.0199;
            obj.rhoS   = 0.2;
            obj.rhoM   = 0.677777777777778;            
            obj.rho1   = 0.9999;
            obj.qRho0  = 32;
            obj.qRhoS  = 11.322033898305085;
            obj.qRhoM  = 3;            
            obj.qRho1  = 25;            
            obj.alpha  = 23.620689655172413;   
            obj.beta   = 1.444444444444444;            
        end
        
        function xi = computeXi(obj)
            xi = SuperEllipseParamsRelator.xi(obj.m1,obj.m2);
        end
        
        function q = initQ(obj)
            nP = length(obj.m1);
            q0 = (obj.qRhoM + obj.qRho0)/2;
            q = q0*ones(nP,1);
        end
        
        function computeExponentAnalytical(obj)
            xi = obj.computeXi();
            q = obj.initQ();            
            error = 1;
            while error > 1e-13
                rho  = obj.computeRho(q);
                qNew = obj.computeQ(xi,rho);
                error = max(abs(q-qNew))/max(abs(q));
                q = qNew;
            end
            obj.value = q;
        end
        
        function rho = computeRho(obj,q)
            rho = SuperEllipseParamsRelator.rho(obj.m1,obj.m2,q);
            rho = max(min(rho,obj.rho1),obj.rho0);
        end
        
        function q = computeQForSymmetricHoles(obj,rho)
            s.rho0  = obj.rho0;
            s.rho1  = obj.rho1;
            s.rhoM  = obj.rhoM;
            s.rhoS  = obj.rhoS;
            s.qRho0 = obj.qRho0;
            s.qRho1 = obj.qRho1;
            s.qRhoM = obj.qRhoM;
            s.qRhoS = obj.qRhoS;
            qSymm = QforSymmetricHolesComputer(s);
            qSymmF = qSymm.compute();
            q = qSymmF(rho);
        end
        
        function q = computeQForElongatedHoles(obj,rho)
            rhoA = (rho - obj.rho0)/(obj.rho1-obj.rho0);
            q = obj.qRho0 + (rhoA).^(obj.beta)*(obj.qRho1-obj.qRho0);
        end
        
        function xiS = computeXiScaled(obj,xi,rho)
            xiMin = obj.computeXiMin(rho);
            xiMax = obj.computeXiMax(rho);
            xiMean = pi/4;
            xiDif  = (xiMax - xiMin)/2;
            xiS = abs(xi - xiMean)./(xiDif);
        end        
        
        function xiMax = computeXiMax(obj,rho)
            qMax = obj.qRho0;            
            xiMax = SuperEllipseParamsRelator.xiFromMxAndRho(obj.m1Max,rho,qMax);
        end
        
        function xiMin = computeXiMin(obj,rho)
            qMax = obj.qRho0;
            xiMin = SuperEllipseParamsRelator.xiFromMyAndRho(obj.m2Max,rho,qMax);         
        end
        
    end
    
end