classdef PostProcessDataBaseCreator < handle

    properties (Access = private)
        mesh 
        outFileName
    end
    
    methods (Access = public)
        
        function obj = PostProcessDataBaseCreator(dI)
            obj.init(dI);
        end

        function data = create(obj)
            d.outFileName = obj.outFileName;
            d.coordinates = obj.mesh.coord;
            d.connectivities = obj.mesh.connec;
            d.nnode = size(obj.mesh.connec,2);
            d.npnod = size(obj.mesh.coord,1);
            d.gtype = obj.mesh.type;
            d.nelem = size(obj.mesh.connec,1);
            d.etype = obj.computeGiDElementType(d.gtype);
            data = d;
        end
        
       
    end

    methods (Access = private)

        function init(obj,cParams)
            if isfield(cParams,'outName')
                cParams.outFileName = cParams.outName;
            end
            obj.mesh        = cParams.mesh;
            obj.outFileName = cParams.outFileName;
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
        

        
    end
    
end