classdef Interpolation_geometry < Interpolation
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
%         order
        geometry_type
    end
    
    methods
        function obj = compute(obj,mesh)
            [obj.isoparametric,obj.order] = Interpolation_geometry.get_isoparametric(mesh);
            obj.xpoints = mesh.coord;
            obj.T = mesh.connec;
            obj.geometry_type = mesh.geometryType;
        end
        
    end
    
    methods (Static)
        function [isoparametric,order] = get_isoparametric(mesh)
            nnode = size(mesh.connec,2);  
            switch mesh.geometryType
                
                    case 'TRIANGLE'

                        switch nnode
                            case 3
                                order = 'LINEAR';
                                isoparametric = Triangle_Linear;
                            case 6
                                order = 'QUADRATIC';
                                isoparametric = Triangle_Quadratic;
                            otherwise
                                error('Invalid nnode for element TRIANGLE.');
                        end
                    case 'Triangle_Linear_Mass'
                        isoparametric=Triangle_Linear_Mass;
                    case 'QUAD'
                        switch nnode
                            case 4
                                order = 'LINEAR';
                                isoparametric = Quadrilateral_Bilinear;
                            case 8
                                order = 'QUADRATIC'
                                isoparametric = Quadrilateral_Serendipity;
                            otherwise
                                error('Invalid nnode for element QUADRILATERAL.');
                        end
                    case 'TETRAHEDRA'
                        isoparametric = Tetrahedra;
                        order = 'LINEAR';
                    case 'HEXAHEDRA'
                        isoparametric = Hexahedra;
                        order = 'LINEAR';
                    otherwise
                        error('Invalid mesh type.')
             end
            
        end
    end
end
