classdef FemTests < testRunner
    
    
    properties (Access = protected)
        FieldOfStudy = 'FEM'
        tests
    end
    
    methods (Access = public)
        function obj = FemTests()
            obj@testRunner();
        end
    end
    
    methods (Access = protected)
        function loadTests(obj)
            obj.tests = {...  
                'test2dStokes_triangle';                                
                'test2dMicro';
                'test2dQuad';                
                'test2dTriangle';
                'test3dTetrahedra';
                'test3dHexahedra';
                };
        end
    end
    
end

