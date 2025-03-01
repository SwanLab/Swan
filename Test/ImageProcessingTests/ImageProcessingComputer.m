classdef ImageProcessingComputer < handle

    properties (Access = public)
        computation
        variables
    end

    properties (Access = private)
        testName
        imageProblem
    end

    methods (Access = public)
        function obj = ImageProcessingComputer(cParams)
            obj.testName = cParams.testName;
        end

        function compute(obj)
            problemName  = obj.testName;
            sLoader      = SettingsLoader(problemName);
            settings     = sLoader.settings;
            imProblem = DenoisingProblem(settings);
            imProblem.solve();
            obj.imageProblem = imProblem;
            obj.variables = imProblem;
            obj.computation = obj;
        end
    end
end