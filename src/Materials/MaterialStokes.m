classdef MaterialStokes < Material

    properties (Access = public)
        mu
        %nElem
    end

    methods (Access = public)

        function obj = MaterialStokes(cParams)
            %obj.nElem = cParams.nelem;
            obj.mu = ones(4,4,cParams.nelem);
            %obj.mu = obj.computeMaterial();
        end
    end

    methods (Access = protected)

        function evaluateNew(obj,~)
        end

    end


    methods (Access = private)

        function m = computeMaterial(obj)
            m = zeros(4,4,obj.nElem);
            m(1,1,:) = 1;
            m(2,2,:) = 1;
            m(3,3,:) = 1;
            m(4,4,:) = 1;
        end

    end

end