classdef MaterialTestComputer < handle

    properties (Access = public)
        comparingMaterialMessage
    end

    properties (Access = private)
        materialOriginal
        materialTest
    end
    
    methods (Access = public)

        function obj = MaterialTestComputer(cParams)
            obj.init(cParams);
        end

        function compute(obj)
            obj.compareMaterial();
            result = obj.comparingMaterialMessage;
            disp(result);
        end

    end


    methods (Access = private)

        function init(obj,cParams)
            obj.saveInput(cParams);
            obj.loadOriginalValues();
        end

        function saveInput(obj,cParams)
            obj.materialTest      = cParams.material;
        end

        function loadOriginalValues(obj)
            load("datas.mat","material");
            obj.materialOriginal   = material;
        end

        function compareMaterial(obj)
           try
                assert(isequal(obj.materialOriginal, obj.materialTest), ...
                    'The material does not have the expected value.');
                obj.comparingMaterialMessage = 'The material has the expected value.';
           catch ME
                obj.comparingMaterialMessage = ME.message;
           end
        end

    end
    
end