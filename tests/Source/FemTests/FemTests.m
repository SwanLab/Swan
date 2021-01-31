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
                'test3dTetrahedra';
                'test3dHexahedra';                
                'test2dMicro';
                'test2dQuad';                
                'test2dTriangle';
                'testPrincipalDirection2D';
                'testPrincipalDirection3D';                
                };
        end
    end
    
end

