classdef FullInnerMeshCreator < handle

    properties (Access = protected) % Inputs
        type
        unfittedMesh
    end
    
    methods (Static, Access = public)

        function ime = create(cParams)
            switch cParams.type
                case 'Matlab'
                    ime = FullInnerMeshCreator_Matlab(cParams);
                case 'GiD'
                    h = cParams.unfittedMesh.innerMesh.mesh.computeMeanCellSize();
                    cParams.meshElementSize = num2str(h);
                    ime = FullInnerMeshCreator_GiD(cParams);
            end
        end
        
    end
    
end