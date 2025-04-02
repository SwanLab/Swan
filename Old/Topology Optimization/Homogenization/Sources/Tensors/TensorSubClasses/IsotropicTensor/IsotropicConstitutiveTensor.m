classdef IsotropicConstitutiveTensor < Stiffness3DTensor
    
    properties (Access = protected)
        E
        nu
        mu
        lambda
        lambda2D
        kappa2D
        kappa
    end
    
    methods (Access = public)
        
        function obj = IsotropicConstitutiveTensor(E,nu)
            obj.createLameParameters(E,nu);
            obj.generateTensor()
        end
        
        function t = clone(obj)
            EV = obj.getYoung();
            nuV = obj.getPoisson();
            t = IsotropicConstitutiveTensor(EV,nuV);
        end
            
        function E = getYoung(obj)
            E = obj.E;
        end
        
        function nu = getPoisson(obj)
            nu = obj.nu;
        end

        function m = getMu(obj)
            m = obj.mu;
        end
        
        function l = getLambda(obj)
            l = obj.lambda;
        end
 
        function l = getLambda2D(obj)
            l = obj.lambda2D;
        end
        
        function k = getKappa(obj)
            k = obj.kappa;
        end

        function k = getKappa2D(obj)
            k = obj.kappa2D;
        end
        
    end
    
    methods (Access = private)
        
        function createLameParameters(obj,E,nu)
            obj.E        = E;
            obj.nu       = nu;
            obj.mu       = obj.computeMuFromYoungAndPoisson(obj.E,obj.nu);
            obj.lambda   = obj.computeLambdaFromYoungAndPoisson(obj.E,obj.nu);
            obj.lambda2D = obj.computeLambda2DFromYoungAndPoisson(obj.E,obj.nu);
            obj.kappa    = obj.computeKappaFromLambdaAndMu(obj.lambda,obj.mu);
            obj.kappa2D  = obj.computeKappa2DFromLambdaAndMu(obj.lambda,obj.mu);
            obj.simplifyExpressions()
        end
        
        function  simplifyExpressions(obj)
            if obj.isDataSymbolic()
               obj.E        = obj.simplifyExpression(obj.E); 
               obj.nu       = obj.simplifyExpression(obj.nu);
               obj.mu       = obj.simplifyExpression(obj.mu);
               obj.lambda   = obj.simplifyExpression(obj.lambda);
               obj.lambda2D = obj.simplifyExpression(obj.lambda2D);
               obj.kappa    = obj.simplifyExpression(obj.kappa);
               obj.kappa2D  = obj.simplifyExpression(obj.kappa2D);
            end            
        end
        
        function itIs = isDataSymbolic(obj)
            isYoungSymbolic   = obj.isSymbolic(obj.E);
            isPoissonSymbolic = obj.isSymbolic(obj.nu);
            itIs = isYoungSymbolic || isPoissonSymbolic;
        end
        
        function a = simplifyExpression(obj,a)
            if obj.isSymbolic(a)
                a = simplify(a);
            end
        end
        
        function generateTensor(obj)
            d = obj.getDimension();
            I = eye(d,d);
            T = zeros(obj.getTensorSize);
            for i = 1:d
                for j = 1:d
                    for k = 1:d
                        for l = 1:d
                            T(i,j,k,l) =  obj.mu*(I(i,k)*I(j,l) + I(i,l)*I(j,k))...
                                + obj.lambda*I(i,j)*I(k,l);
                        end
                    end
                end
            end
            obj.setValue(T);
        end
        
    end
    
    methods (Access = private, Static)
        
        function itIs = isSymbolic(a)
            itIs = isa(a,'sym');
        end
    end
    
    methods (Access = public, Static)
        
        function  tens = createWithLambdaAndMu(lambda,mu)
            E  = IsotropicConstitutiveTensor.computeYoungFromLambdaAndMu(lambda,mu);
            nu = IsotropicConstitutiveTensor.computePoissonFromLambdaAndMu(lambda,mu);
            tens = IsotropicConstitutiveTensor(E,nu);
        end
        
        function E = computeYoungFromLambdaAndMu(lambda,mu)
            E = mu*(3*lambda+2*mu)/(lambda+mu);
        end
        
        function nu = computePoissonFromLambdaAndMu(lambda,mu)
            nu = lambda/(2*(lambda + mu));
        end
        
        function lambda2D = computeLambda2DFromYoungAndPoisson(E,nu)
            lambda2D = E*nu/(1+nu)/(1-nu);
        end
        
        function mu = computeMuFromYoungAndPoisson(E,nu)
            mu = (E/(2*(1+nu)));
        end
        
        function lambda = computeLambdaFromYoungAndPoisson(E,nu)
            lambda = (E*nu/((1+nu)*(1-2*nu)));
        end
        
        function kappa2D = computeKappa2DFromLambdaAndMu(lambda,mu)
            kappa2D = lambda + mu;
        end
        
        function kappa = computeKappaFromLambdaAndMu(lambda,mu)
            kappa = lambda + 2/3*mu;
        end
        
    end
    
end

