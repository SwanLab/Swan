classdef MeshCreatorFromRVE < handle

    methods (Access = public, Static)

        function obj = create(cParams)
            switch cParams.meshReference.geometryType
                case 'Surface'
                   obj = MeshCreatorFromRVE2D(cParams);
                case 'Volume'
                   obj = MeshCreatorFromRVE3D(cParams);
            end
        end

    end

end