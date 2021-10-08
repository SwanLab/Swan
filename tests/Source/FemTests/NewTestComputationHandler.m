classdef (Abstract) NewTestComputationHandler < handle
    
    methods (Static)

        function computer = create(comp_type)    
            switch comp_type
                case {'FEM_SOLVER'}
                    computer = NewFemSolver();
                case {'ITERATIVE'}
                    computer = IterativeSolver();
                otherwise
                    error('Invalid Computer Type.')
            end
        end

    end
    
    methods (Abstract, Access = public)
        compute
    end

end