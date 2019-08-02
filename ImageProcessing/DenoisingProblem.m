classdef DenoisingProblem < handle
    
    properties (GetAccess = public, SetAccess = private)
        optimizedImage
    end
    
    properties (Access = private)
        originalImage
        imageSize
        discreteGradient
        noisyImage
        noisyImageNorm
        
        
        En
        gap
        maxIter
        lipschitzConstant
        lambda
        pDual
    end
    
    methods (Access = public)
        
        function obj = DenoisingProblem(cParams)
            obj.init(cParams);
            obj.createDiscreteGradient();
            obj.createNoisyImage(cParams);
            obj.computeNoisyImageNorm();
        end
        
        function solve(obj)
            mn = obj.imageSize.rowsTimesColumns;     
            obj.pDual = zeros(2*mn,1);
            for i=1:obj.maxIter
                obj.computeGradientStep();
                obj.projectInTheBall();
                obj.computeU();
                obj.En(i) = obj.computeEnergy();
                obj.gap(i) = obj.computeDualGap();
            end
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            im = obj.readImage(cParams);
            obj.computeImageSize(im);
            obj.transformImageInVectorForm(im)
            obj.maxIter = cParams.maxIter;
            obj.lipschitzConstant = cParams.lipschitzConstant;
            obj.lambda = cParams.totalVariationWeigth;
        end
        
        function computeImageSize(obj,image)
            [m,n] = size(image);
            obj.imageSize.rows = m;
            obj.imageSize.columns = n;
            obj.imageSize.rowsTimesColumns = m*n;
        end
        
        function transformImageInVectorForm(obj,image)
            obj.originalImage = image(:);
        end
        
        function createDiscreteGradient(obj)
            m  = obj.imageSize.rows;
            n  = obj.imageSize.columns;
            mn = obj.imageSize.rowsTimesColumns;
            I  = reshape(1:m*n,m,n);
            east  = [I(:,2:end), I(:,end)];
            north = [I(2:end,:); I(end,:)];
            D1 = sparse(I,east,1,mn,mn)  -speye(mn,mn);
            D2 = sparse(I,north,1,mn,mn) -speye(mn,mn);
            obj.discreteGradient = [D1 ; D2];
        end
        
        function createNoisyImage(obj,cParams)
            u0    = obj.originalImage;
            L     = cParams.lipschitzConstant;
            a     = cParams.noiseAmplitud;
            sigma = a/L;
            uN    = u0 + sigma*randn(size(u0));
            obj.noisyImage = uN;
        end
        
        function computeNoisyImageNorm(obj)
            g = obj.noisyImage;
            obj.noisyImageNorm = 0.5*(g'*g);
        end
        
        function computeGradientStep(obj)
            L = obj.lipschitzConstant;
            tauV = 1/L;
            D = obj.discreteGradient;
            g = obj.noisyImage;
            p = obj.pDual;
            p = p - tauV*D*(D'*p - g);
            obj.pDual = p;
        end
        
        function projectInTheBall(obj)
            lam = obj.lambda;
            normP = obj.computeNormP();
            no  = max(1,normP/lam);
            p = obj.pDual;
            p = p./[no;no];
            obj.pDual = p;
        end
        
        function normP = computeNormP(obj)
            p  = obj.pDual;
            mn = obj.imageSize.rowsTimesColumns;
            normP = hypot(p(1:mn),p(mn+1:end));
        end
        
        function computeU(obj)
            g = obj.noisyImage;
            D = obj.discreteGradient;
            p = obj.pDual;
            obj.optimizedImage = g(:) - D'*p;
        end
        
        function E = computeEnergy(obj)
            p = obj.pDual;
            D = obj.discreteGradient;
            g = obj.noisyImage;
            gN = obj.noisyImageNorm;
            r = D'*p - g;
            E = 0.5*(r'*r)/gN;
        end
        
        function gap = computeDualGap(obj)
            lam = obj.lambda;
            ut  = obj.optimizedImage;
            p   = obj.pDual;
            D   = obj.discreteGradient;
            gN  = obj.noisyImageNorm;
            q   = D'*p;
            gap = lam*sum(abs(D*ut)) + q'*ut;
            gap = gap/gN;
        end
        
    end
    
    methods (Access = private, Static)
        
        function im = readImage(cParams)
            image = cParams.imageFile;
            im = double(imread(image));
        end
        
    end
    
end