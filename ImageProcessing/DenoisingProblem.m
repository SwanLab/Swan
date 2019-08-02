classdef DenoisingProblem < handle
    
    properties (Access = private)
        originalImage
        imageSize
        discreteGradient
        noisyImage
        noisyImageNorm
        optimizedImage
    end
    
    methods (Access = public)
        
        function obj = DenoisingProblem(cParams)
            obj.init(cParams);
            obj.createDiscreteGradient();
            obj.createNoisyImage(cParams);
        end
        
        function solve(obj)
            
        end
        
    end
    
    
    methods (Access = private)
        
        function init(obj,cParams)
            im = obj.readImage(cParams);
            obj.computeImageSize(im);
            obj.transformImageInVectorForm(im)
        end
        
        function computeImageSize(obj,image)
            [m,n] = size(image);
            obj.imageSize.rows = m;
            obj.imageSize.columns = n;
            obj.imageSize.rowsTimesColums = m*n;
        end
        
        function transformImageInVectorForm(obj,image)
            obj.originalImage = image(:);
        end
        
        function createDiscreteGradient(obj)
            m  = obj.imageSize.rows;
            n  = obj.imageSize.columns;
            mn = obj.imageSize.rowsTimesColums;
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
        
    end
    
    methods (Access = private, Static)
        
        function im = readImage(cParams)
            image = cParams.imageFile;
            im = double(imread(image));
        end
        
    end
    
end