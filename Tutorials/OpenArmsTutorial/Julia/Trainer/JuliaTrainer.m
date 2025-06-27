classdef JuliaTrainer < handle
    properties (Access = public)
        id
        type
        objectiveFunction
        designVariable
    end

    properties (Access = private)
        data  % Struct returned from Julia Create() method (e.g. might include metadata)
    end
    
    methods (Access = public, Static)
        
     function obj = create(params)
            % Call the Julia Create function in Trainer module
            result = callJuliaClass('Trainer', 'Create', params);

            % Create object and assign metadata
            obj = JuliaTrainer();  % No-arg constructor
            obj.id = result.id;
            obj.type = result.type;

            if isfield(params, 'costFunc')
                obj.objectiveFunction = params.costFunc;
            end
            if isfield(params, 'designVariable')
                obj.designVariable = params.designVariable;
            end
        end

    end
end
