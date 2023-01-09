classdef STLExporter < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private) % Inputs
        filename
        type
        mesh
    end
    
    properties (Access = private) % Calculated
        triMesh
    end
    
    methods (Access = public)
        
        function obj = STLExporter(cParams)
            obj.init(cParams)
        end

        function export(obj)
            file = obj.filename;
            switch obj.type
                case 'Matlab'
                    obj.exportUsingMatlab(file);
                case 'GiD'
                    obj.exportUsingGiD(file);
            end
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.mesh     = cParams.mesh;
            obj.type     = cParams.type;
            obj.filename = cParams.filename;
        end

        function triangulateMesh(obj)
            P = obj.mesh.coord;
            T = obj.mesh.connec;
            obj.triMesh = triangulation(T,P);
        end

        function exportUsingMatlab(obj, file)
            obj.triangulateMesh();
            stlwrite(obj.triMesh, [file '.stl'])
        end

        function exportUsingGiD(obj)
        end
        
    end
    
end