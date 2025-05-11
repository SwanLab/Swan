classdef ImageProcessingComputer < handle

    properties (Access = public)
        computation
        variables
    end

    properties (Access = private)
        imageFile
        lipschitzConstant
        totalVariationWeight
        noiseAmplitud
        maxIter
        optimizer
    end

    properties (Access = private)
        imageProblem
    end

    methods (Access = public)
        function obj = ImageProcessingComputer(cParams)
            obj.imageFile            = cParams.imageFile;
            obj.lipschitzConstant    = cParams.lipschitzConstant;
            obj.totalVariationWeight = cParams.totalVariationWeight;
            obj.noiseAmplitud        = cParams.noiseAmplitud;
            obj.maxIter              = cParams.maxIter;
            obj.optimizer            = cParams.optimizer;
        end

        function compute(obj)
            s.imageFile            = obj.imageFile;
            s.totalVariationWeight = obj.totalVariationWeight;
            s.lipschitzConstant    = obj.lipschitzConstant;
            s.noiseAmplitud        = obj.noiseAmplitud;
            s.maxIter              = obj.maxIter;
            s.optimizer            = obj.optimizer;
            imProblem              = DenoisingProblem(s);
            imProblem.solve();
            obj.imageProblem = imProblem;
            obj.variables    = imProblem;
            obj.computation  = obj;
        end
    end
end