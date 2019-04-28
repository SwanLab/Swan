classdef TopOptTests < testRunner

    properties (Access = protected)
        FieldOfStudy = 'Topology Optimization'
        tests
    end

    methods (Access = public)

        function obj = TopOptTests()
            obj@testRunner();
        end

    end

    methods (Access = protected)

        function loadTests(obj)
            obj.tests = {... 
            'testProjectedSlerp';                                                                   
            'testCantilever3';                
            'testBridge2';                                      
            'testMicro';                
            'testCantilever';   
            'testCantilever2';                        
            'testGripping';                
            'testBridge';           
            'testMicro2';            
            };
        end

    end

end
