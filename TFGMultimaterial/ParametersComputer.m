classdef ParametersComputer < handle

    properties (Access = public)
        params
    end

    methods (Access = public)

        function obj = ParametersComputer()
            obj.init();
        end

    end

    methods (Access = private)

        function init(obj)
            load viga2x1.mat;
            obj.params.kmin = eps; % minimum allowed 'k'. Used by the line-search procedure
            obj.params.stop = 1.0*pi/180; % stop criterion
            obj.params.penalization = 3; % Augmented Lagrangian penalization
            obj.params.penalty = [0.0 0.0 0.0];   % method 1,2 and 3
            obj.params.volfrac = 0.6*[1/3 1/3 1/3];   % method 2 and 3
            obj.params.voleps  = [0.01 0.01 0.01];  % method 2 and 3
            obj.params.auglag  = [0.05 0.05 0.05];   % method 3
            obj.params.epsilon = 0.01;  % method 3       
        end
    end

end