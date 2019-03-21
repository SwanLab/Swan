classdef ExploringSettingsTests < testRunner
    
    properties (Access = protected)
        FieldOfStudy = 'ExploringSettingsTests'
        tests
    end
    
    methods (Access = public)
        
        function obj = ExploringSettingsTests()
            obj@testRunner();
        end
        
    end
    
    methods (Access = protected)
        
        function loadTests(obj)
            obj.tests = {...
                'testExplorSettNumHomogDefault';
                'testExplorSettNumHomogCustom';
                'testExplorSettLevSetDefault';
                'testExplorSettLevSetCustom';
                };
            
        end
        
    end
    
end

