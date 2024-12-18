classdef ProjectorFactory < handle

    methods (Static, Access = public)

        function obj = create(cParams)
            switch cParams.projectorType
                case 'P0'
                    obj = Projector_toP0(cParams);
                case {'P1','P2','P3','P1D'}
                    obj = Projector_toLagrangian(cParams);
                case 'RT'
                    obj = Projector_toRaviartThomas(cParams);
                case 'N'
                    obj = Projector_toNedelec(cParams);
                case 'RigidBody'
                    obj = Projector_toRigidBody(cParams);
                case 'ModalFunction'
                    obj = Projector_toModalFunction(cParams);
            end
        end

    end
end
