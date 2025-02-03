classdef StokesProblemSolverTestComputer < handle

    properties (Access = public)
        comparingVelocityFunMessage
        comparingPressureFunMessage
    end

    properties (Access = private)
        
    end

    properties (Access = private)
        velocityFunOriginal
        pressureFunOriginal
        velocityFunTest
        pressureFunTest
    end
    
    methods (Access = public)

        function obj = StokesProblemSolverTestComputer(cParams)
            obj.init(cParams);
        end

        function compute(obj)
            obj.compareVelocityFun();
            obj.comparePressureFun();
            result = [obj.comparingVelocityFunMessage, newline, obj.comparingPressureFunMessage];
            disp(result);
        end

    end


    methods (Access = private)

        function init(obj,cParams)
            obj.saveInput(cParams);
            obj.loadOriginalValues();
        end

        function saveInput(obj,cParams)
            obj.velocityFunTest          = cParams.velocityFun;
            obj.pressureFunTest          = cParams.pressureFun;
        end

        function loadOriginalValues(obj)
            load("datas.mat","velocityFun","pressureFun");
            obj.velocityFunOriginal   = velocityFun;
            obj.pressureFunOriginal   = pressureFun;
        end

        function compareVelocityFun(obj)
           try
                assert(isequal(obj.velocityFunOriginal, obj.velocityFunTest), ...
                    'The velocity function does not have the expected value.');
                obj.comparingVelocityFunMessage = 'The velocity function has the expected value.';
           catch ME
                obj.comparingVelocityFunMessage = ME.message;
           end
        end

        function comparePressureFun(obj)
           try
                assert(isequal(obj.pressureFunOriginal, obj.pressureFunTest), ...
                    'The pressure function does not have the expected value.');
                obj.comparingPressureFunMessage = 'The pressure function has the expected value.';
           catch ME
                obj.comparingPressureFunMessage = ME.message;
           end
        end

    end
    
end