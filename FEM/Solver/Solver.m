classdef Solver < handle

    methods (Static)

        function stype = create(cParams)
            %             solver_type = 'DIRECT';
            switch cParams.solverType
                case {'DIRECT'}
                    stype = Direct_Solver();

                case {'CHOLEVSKY'}
                    % At least up to ndof ~5e4, Direct is still faster
                    stype = Cholesky_Direct_Solver();

                case {'ITERATIVE'}
%                     error('Not implemented yet')
                    switch cParams.iterativeSolverType
                        case{'JACOBI'}
                            stype = JacobiSolver(cParams);
                        case{'GAUSS'}
                            stype = GaussSolver(cParams);
                        case{'PCG'}
                            stype = PCG(cParams);
                        case{'CG'}    
                            stype = ConjugateGradientSolver(cParams);    
                        case{'MULTIGRID'}
                            stype = Multigrid(cParams);
                    end

%                 case 'PCG'
%                     % Preconditioned conjugate gradient
%                     stype = PCG(cParams);

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