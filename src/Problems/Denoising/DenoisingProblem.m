classdef DenoisingProblem < handle
    
    properties (GetAccess = public, SetAccess = private)
        optimizedImage
        plottingData
    end
    
    properties (Access = private)
        originalImage
        imageSize
        discreteGradient
        noisyImage
        noisyImageNorm
        
        optimizer
        designVariable
        iterator
        
        
        
        energyFunction
        l1Proximal
        dualGap
        
        lambda
    end
    
    methods (Access = public)
        
        function obj = DenoisingProblem(cParams)
            obj.init(cParams);
            obj.createDiscreteGradient();
            obj.createNoisyImage(cParams);
            obj.computeNoisyImageNorm();
            obj.createDesignVariable();
            obj.createIterator(cParams);
            obj.createEnergyFunction(cParams);
            obj.createL1Proximal(cParams);
            obj.createOptimizer(cParams);
        end
        
        function solve(obj)
            while obj.iterator.hasNotFinished()
                obj.optimizer.update();
                obj.designVariable.update();
                obj.energyFunction.computeCost();
                obj.computeU();
                obj.computeDualGap();
                obj.updatePlottingData();
                obj.iterator.update;
            end
        end
        
    end
    
    methods (Access = private)
        
        function updatePlottingData(obj)
            p = obj.plottingData;
            i = obj.iterator.value;
            p.cost(i) = obj.energyFunction.value;
            p.dualGap(i) = obj.dualGap; 
            obj.plottingData = p;
        end
        
        function init(obj,cParams)
            im = obj.readImage(cParams);
            obj.computeImageSize(im);
            obj.transformImageInVectorForm(im)
            obj.lambda = cParams.totalVariationWeight;
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
        
        function createDesignVariable(obj)
            mn = 2*obj.imageSize.rowsTimesColumns;
            s.xLength = mn;
            obj.designVariable = DesignImagVariable(s);
        end
    
        function computeU(obj)
            g = obj.noisyImage;
            D = obj.discreteGradient;
            p = obj.designVariable.value;
            obj.optimizedImage = g(:) - D'*p;
        end

        function createEnergyFunction(obj,cParams)
            s.lipschitzConstant = cParams.lipschitzConstant;
            s.A = obj.discreteGradient';
            s.b = obj.noisyImage;
            s.designVariable = obj.designVariable;
            c = QuadraticFunction(s);
            obj.energyFunction = c;
        end
        
        function createL1Proximal(obj,cParams)
            s.lambda = cParams.totalVariationWeight;
            s.imageSize = obj.imageSize;
            s.designVariable = obj.designVariable;
            obj.l1Proximal = L1VectorNormProximal(s);
        end
        
        function createIterator(obj,cParams)
            s.maxIter = cParams.maxIter;
            obj.iterator = Iterator(s);
        end
        
        function computeDualGap(obj)
            lam = obj.lambda;
            ut  = obj.optimizedImage;
            p   = obj.designVariable.value;
            D   = obj.discreteGradient;
            gN  = obj.noisyImageNorm;
            q   = D'*p;
            gap = lam*sum(abs(D*ut)) + q'*ut;
            obj.dualGap = gap/gN;
        end
        
        function createOptimizer(obj,cParams)
            sg.designVariable = obj.designVariable;
            sg.differentiableFunction = obj.energyFunction;
            sg.designVariable = obj.designVariable;
            s.gradientMethodParams = sg;
            
            sm.designVariable = obj.designVariable;
            sm.iterator       = obj.iterator;
            s.momentumParams  = sm;
            
            s.proximal = obj.l1Proximal;
            s.type     = cParams.optimizer;
            s.iterator = obj.iterator;
            obj.optimizer = SplittingAlgorithm.create(s);
        end
        
    end
    
    methods (Access = private, Static)
        
        function im = readImage(cParams)
            image = cParams.imageFile;
            im = double(imread(image));
        end
        
    end
    
end