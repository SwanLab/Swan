classdef ReadingFilesTests < testRunner
    
    
    properties (Access = protected)
        FieldOfStudy = 'ReadingFiles'
        tests
    end
    
    methods (Access = public)
        function obj = ReadingFilesTests()
            obj@testRunner();
        end
    end
    
    methods (Access = protected)
        function loadTests(obj)
            obj.tests = {...  
                'testReadingGmsh';                
                };
        end
    end
    
end