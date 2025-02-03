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
            obj.tolerance = 1e-6;
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
                assert(all(abs(obj.LiftOriginal - obj.LiftTest) < obj.tolerance, 'all'), ...
                'The Lift forces do not have the expected value.');
                 obj.comparingLiftMessage = 'The Lift forces have the expected value.';
            catch ME
                obj.comparingLiftMessage = ME.message;
            end
        end

        function compareDrag(obj)
            try
                assert(all(abs(obj.DragOriginal - obj.DragTest) < obj.tolerance, 'all'), ...
                    'The Drag forces do not have the expected value.');
                obj.comparingDragMessage = 'The Drag forces have the expected value.';
            catch ME
                obj.comparingDragMessage = ME.message;
            end
        end

    end
    
end