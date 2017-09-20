classdef Preprocess<handle
    % Only reads
    
    properties
    end
    
    methods(Static)
        function [coord, connec] = readFromGiD()
            coord = [
                1            0            0            0
                2         0.01         0.01            0
                3            0         0.02            0
                4         0.02            0            0
                5         0.02         0.02            0
                6         0.01         0.03            0
                7         0.03         0.01            0
                8            0         0.04            0
                9         0.04            0            0
                10        0.03         0.03            0
                11        0.04         0.02            0
                12        0.02         0.04            0
                13        0.04         0.04            0
                ];
            
            %%
            connec = [
                1   5   2   4   0
                2   3   2   5   0
                3   1   2   3   0
                4   4   2   1   0
                5   12  6   5   0
                6   8   6   12  0
                7   3   6   8   0
                8   5   6   3   0
                9   11  7   9   0
                10  5   7   11  0
                11  4   7   5   0
                12  9   7   4   0
                13  13  10  11  0
                14  12  10  13  0
                15  5   10  12  0
                16  11  10  5   0
                ];
            
        end
        function [fixnodes, forces] = getBC()
            % Dirichlet
            % Node - Dimension - Value
            fixnodes = [
                1 1 0
                1 2 0
                3 1 0
                3 2 0
                8 1 0
                8 2 0
                ];
            
            % Neumann --> Fpunc (global)
            forces = [
                11 2 -1
                ];
        end
    end
    
end


