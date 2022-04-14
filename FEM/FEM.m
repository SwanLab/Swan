classdef FEM < handle
    
    properties (Access = protected)
        inputReader
    end
    
    methods (Access = public)
        
        function obj = FEM()
            obj.inputReader = FemInputReader_GiD();
        end
        
    end
    
    methods (Static, Access = public)
        
        function obj = create(s)
            switch s.type
                case 'ELASTIC'
                    switch s.scale
                        case 'MACRO'
                            obj = ElasticProblem(s);
                        case 'MICRO'
                            obj = ElasticProblemMicro(s);
                    end
                case 'THERMAL'
                    obj = ThermalProblem(s);
                case 'DIFF-REACT'
                    obj = DiffReactProblem.create(s);
                case 'HYPERELASTIC'
                    obj = Hyperelastic_Problem(s);
                case 'Stokes'
                    obj = Stokes_Problem(fileName);
            end
        end
        
    end

end
