classdef (Abstract) TestComputer < handle
    
    properties (Abstract)
        computation
    end
    
    methods (Static)

        function computer = create(solver_type, s)    
            switch solver_type
                case {'STOKES'}
                    computer = StokesComputer(s);
                case {'FEM'}
                    computer = FemComputer(s);
                case {'MICRO'}
                    computer = MicroComputer(s);
                case {'TOPOPT'}
                    computer = TopOptComputer(s);
                otherwise
                    error('Invalid Computer Type.')
            end
        end

    end
    
    methods (Abstract, Access = public)
        compute
    end

end