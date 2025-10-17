classdef Material < BaseFunction

    methods (Access = public, Static)

        function material = create(cParams)
            f = MaterialFactory();
            material = f.create(cParams);
        end

    end

    methods (Access = protected)

        function init(obj,cParams)
            obj.mesh  = cParams.mesh;
            obj.ndimf = obj.mesh.ndim.*ones(1,4);
            obj.ndimfTotal = prod(obj.ndimf);
        end

    end

end