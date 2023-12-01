classdef Preconditioner < handle

    methods (Static)

        function ptype = create(cParams)

            switch cParams.preconditionerType
                case {'JACOBI'}
                    ptype = JacobiPreconditioner(cParams);

                case {'CHOLEVSKY'}
                    % At least up to ndof ~5e4, Direct is still faster
                    ptype = Cholesky_Direct_Solver();

                otherwise
%                     error('Invalid preconditioner type.')
            end
        end

    end

end