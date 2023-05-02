classdef LHSFactory < LHSComputer

    methods (Access = public, Static)

        function obj = create(cParams)
            switch cParams.funcType
                case 'Symmetric'
                    obj = SymmetricLHSComputer(cParams);
                case 'Unsymmetric'
                    obj = UnsymmetricLHSComputer(cParams);
            end
        end
    end

end