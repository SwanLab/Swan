classdef AeroForcesCalculationTestComputer < handle

    properties (Access = public)
        comparingLiftMessage
        comparingDragMessage
    end

    properties (Access = private)
        tolerance
    end

    properties (Access = private)
        LiftOriginal
        DragOriginal
        LiftTest
        DragTest
    end
    
    methods (Access = public)

        function obj = AeroForcesCalculationTestComputer(cParams)
            obj.init(cParams);
        end

        function compute(obj)
            obj.compareLift();
            obj.compareDrag();
            result = [obj.comparingLiftMessage, newline, obj.comparingDragMessage];
            disp(result);
        end

    end


    methods (Access = private)

        function init(obj,cParams)
            obj.saveInput(cParams);
            obj.loadOriginalfinalMesh();
            obj.tolerance = 1e-4;
        end

        function saveInput(obj,cParams)
            obj.LiftTest    = cParams.L;
            obj.DragTest    = cParams.D;
        end

        function loadOriginalfinalMesh(obj)
            load("datas.mat","Li","Di");
            obj.LiftOriginal = Li;
            obj.DragOriginal = Di;
        end

        function compareLift(obj)
            try
                relativeError = abs(obj.LiftOriginal - obj.LiftTest) / abs(obj.LiftOriginal);
                assert(relativeError < obj.tolerance, ...
                    sprintf('The lift force does not match the expected value.\nRelative error: %.4f', relativeError));
                obj.comparingLiftMessage = 'The lift force corresponds to the expected value.';
            catch ME
                obj.comparingLiftMessage = ME.message;
            end
        end

        function compareDrag(obj)
            try
                relativeError = abs(obj.DragOriginal - obj.DragTest) / abs(obj.DragOriginal);
                assert(relativeError < obj.tolerance, ...
                    sprintf('The drag force does not match the expected value.\nRelative error: %.4f', relativeError));
                obj.comparingDragMessage = 'The drag force corresponds to the expected value.';
            catch ME
                obj.comparingDragMessage = ME.message;
            end
        end

    end
    
end