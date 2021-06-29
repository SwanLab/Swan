classdef PostProcessDataBaseCreator < handle
    
    properties (Access = protected)
        data
    end
    
    methods (Access = public)
        
        function obj = PostProcessDataBaseCreator(dI)
            mesh          = dI.mesh;
            d.outFileName = dI.outName;
            d.coordinates = mesh.coord;
            d.connectivities = mesh.connec;
            d.nnode = size(mesh.connec,2);
            d.npnod = size(mesh.coord,1);
            d.gtype = mesh.type;
            d.pdim  = dI.pdim;
            d.nelem = size(mesh.connec,1);
            d.ptype = dI.ptype;
            d.ndim  = obj.computeNdim(d.pdim);
            d.etype = obj.computeGiDElementType(d.gtype);
            obj.data = d;
        end
                
        function d = getValue(obj)
            d = obj.data;
        end
        
    end
       
    methods (Access = private, Static)
        
        function et = computeGiDElementType(gt)
            switch  gt
                case 'TRIANGLE'
                    et = 'Triangle';
                case 'QUAD'
                    et = 'Quadrilateral';
                case 'TETRAHEDRA'
                    et = 'Tetrahedra';
                case 'HEXAHEDRA'
                    et = 'Hexahedra';
            end
        end
        
        function n = computeNdim(p)
            switch p
                case '2D'
                    n = 2;
                case '3D'
                    n = 3;
            end
        end
        
    end
    
end