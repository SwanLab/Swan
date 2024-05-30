classdef VolumeComputer < handle

    properties (Access = private)
        area
        tfi
    end

    methods (Access = public)

        function obj = VolumeComputer(cParams)
            obj.init(cParams);
        end

        function vol = computeVolume(obj)
            surf = obj.area;
            charfunc = obj.tfi;

            vol = surf*charfunc'; %calculate volume assigned to each material including void -- mixed formulation aproach
          % vol  = area*tchi';  %calculate volume assigned to each material including void -- P1 projection aproach
        end
    end

    methods (Access = public)
        
        function init(obj,cParams)
            %obj.area = cParams.mesh.area;
            obj.area = cParams.area;
            obj.tfi = cParams.tfi; 
        end
    end
end