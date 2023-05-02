classdef RHSFactory < RHSComputer

    methods (Access = public, Static)
        function obj = create(cParams)
            switch cParams.funcType
                case 'Symmetric'
                    obj = SymmetricRHSComputer(cParams);
                case 'Unsymmetric'
                    obj = UnsymmetricRHSComputer(cParams);
            end
        end
    end
end