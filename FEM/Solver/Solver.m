classdef Solver < handle

    methods (Static)

        function stype = create(cParams)
%             solver_type = 'DIRECT';
            switch cParams.type
                case {'DIRECT'}
                    stype = Direct_Solver();

                case {'CHOLEVSKY'}
                    % At least up to ndof ~5e4, Direct is still faster
                    stype = Cholesky_Direct_Solver();

                case {'ITERATIVE'}
                    error('Not implemented yet')

                case 'CG'
                    stype = CGsolver();

                case 'Nonlinear'
                    stype = NonLinear_Solver(cParams);

                otherwise
                    error('Invalid solver type.')
            end
        end

    end

end