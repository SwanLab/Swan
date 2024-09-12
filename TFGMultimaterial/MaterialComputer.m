classdef MaterialComputer < handle

    properties (Access = public)
        material
    end

    properties (Access = private)
        charFunc
        designVariable
        p
        t
        m
        materialInterpolator
    end

    methods (Access = public)

        function obj = MaterialComputer(cParams)
            obj.init(cParams);
            %obj.getCharacteristicFunction();
            obj.createMaterial();
        end
    end

    methods (Access = private)

        function init(obj,cParams)
            obj.designVariable = cParams.designVariable;
            obj.p = cParams.p;
            obj.t = cParams.t;
            obj.m = cParams.mesh;
            obj.materialInterpolator = cParams.materialInterpolator;
            obj.charFunc = cParams.charFunc;
        end

       % function getCharacteristicFunction(obj)
            % s.m = obj.m;
            % s.t = obj.t;
            % s.p = obj.p;
            % s.designVariable = obj.designVariable;
            % 
            % charfunc = CharacteristicFunctionComputer(s);
            % [~,tfi] = charfunc.computeFiandTfi();
            % obj.charFunc = tfi;
        %end

        function createMaterial(obj)
            s.type                 = 'DensityBased';
            s.density              = obj.charFunc;
            s.materialInterpolator = obj.materialInterpolator;
            s.dim                  = '2D';
            obj.material = Material.create(s);
        end
    end
end