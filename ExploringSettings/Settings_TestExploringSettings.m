classdef Settings_TestExploringSettings < SettingsNumericalHomogenizer
    
    %  ____                      
    % / ___|_      ____ _ _ __   
    % \___ \ \ /\ / / _` | '_ \  
    %  ___) \ V  V / (_| | | | | 
    % |____/ \_/\_/ \__,_|_| |_|
    %               Topology Optimization Laboratory
    %
    % Specific implementation of NumHomog settings for a specific case
    % Here we set the specific parameters that configure our test
    % This kind of files would replace our current Benchmark files
    
    properties
    end
    methods
        function obj = Settings_TestExploringSettings()
            obj.outFileName = 'test2d_micro';
            obj.testName = 'test2d_micro.m';
            obj.print = true;
            obj.volumeShFuncParams.filterParams.optimizer = 'MMA';
            obj.volumeShFuncParams.filename = 'test2d_micro.m';
            obj.volumeShFuncParams.ptype = 'MICRO';
        end
    end
end