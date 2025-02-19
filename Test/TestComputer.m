classdef (Abstract) TestComputer < handle
    
    properties (Abstract)
        computation
    end
    
    methods (Static)

        function computer = create(solver_type, s)
            switch solver_type
                case {'FEM'}
                    computer = FemComputer(s);
                case {'GMSH'}
                    computer = GMSHComputer(s);
                case {'IMAGE'}
                    computer = ImageProcessingComputer(s);
                case {'MICRO'}
                    computer = MicroComputer(s);
                case {'STOKES'}
                    computer = StokesComputer(s); % Li passem a la classe StokesComputer s (cParams), però com no li demanem res més, només fa el primer mètode
                case {'THERMAL'}
                    computer = ThermalComputer(s);
                case {'TOPOPT'}
                    computer = TopOptComputer(s);
                case {'PROJECT'}
                    computer = ProjectorComputer(s);
                case {'ACADEMIC'}
                    computer = AcademicProblem(s);
                otherwise
                    error('Invalid Computer Type.')
            end
        end

    end
    
    methods (Abstract, Access = public)
        compute
    end

end