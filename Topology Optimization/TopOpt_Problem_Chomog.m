classdef TopOpt_Problem_Chomog < TopOpt_Problem
    properties
    end
    methods
        function obj=TopOpt_Problem_Chomog(settings)
            obj@TopOpt_Problem(settings);
            switch settings.ptype
                case 'Chomog_alphabeta'
                    obj.cost=ShFunc_Chomog_alphabeta(settings);
                otherwise
                    fprintf('Micro problem not implemented');
            end            
            obj.constraint=ShFunc_Volume(settings);
        end
    end
end
