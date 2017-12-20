classdef PhysicalVariables
    %PhysicalVariables Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods (Access = ?Physical_Problem, Static)
        function physicalVars = create(ptype,pdim)
            switch ptype
                case 'ELASTIC'
                    switch pdim
                        case '2D'
                            physicalVars = PhysicalVars_Elastic_2D;
                        case '3D'
                            physicalVars = PhysicalVars_Elastic_3D;
                    end
                case 'THERMAL'
                    error('Still not implemented.')
                otherwise
                    error('Invalid ptype.')
            end
        end
    end
    
end

