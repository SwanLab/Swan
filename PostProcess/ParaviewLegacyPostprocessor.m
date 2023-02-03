classdef ParaviewLegacyPostprocessor < handle
    
    % Its output is a legacy .vtk file. It works for simple fields, but it
    % cannot represent data at Gaussian points, as it requires the newer
    % XML implementation.
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        mesh
        coord
        connec
        filename
        outputFile
        d_u
    end
    
    methods (Access = public)
        
        function obj = ParaviewLegacyPostprocessor(cParams)
            obj.init(cParams);
            obj.openFile();
            obj.writeHeader();
            obj.writeCoords();
            obj.writeConnec();
            obj.writeData();
            obj.saveFile();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj, cParams)
            obj.mesh     = cParams.mesh;
            obj.coord    = cParams.mesh.coord;
            obj.connec   = cParams.mesh.connec;
            obj.filename = cParams.filename;
            obj.d_u = cParams.d_u;
        end
        
        function openFile(obj)
            fullfile = strcat(obj.filename, '.vtk');
            obj.outputFile = fopen(fullfile,'w');
        end
        
        function writeHeader(obj)
            header = '# vtk DataFile Version 2.0 \nSample test \nASCII\n';
            fprintf(obj.outputFile, header);
        end
        
        function writeCoords(obj)
            nnodes = size(obj.coord,1);
            head = '\nDATASET UNSTRUCTURED_GRID\nPOINTS %d float\n';
            fprintf(obj.outputFile, head, nnodes);
            if (size(obj.coord,2) == 2)
                obj.coord = [obj.coord, zeros(nnodes,1)];
            end
            fprintf(obj.outputFile,'%.4f %.4f %.4f\n', obj.coord');
        end
        
        function writeConnec(obj)
            nelems = size(obj.connec,1);
            nnodEl = size(obj.connec,2);
            ndata = nelems * (nnodEl+1);
            headCells = '\nCELLS %d %d\n';
            fprintf(obj.outputFile, headCells, nelems, ndata);
            connecP = [nnodEl*ones(nelems,1), obj.connec-1];
            format = ['%d', repmat(' %d',1,nnodEl),'\n'];
            fprintf(obj.outputFile,format, connecP'); %d fixed!!
            type = obj.getCellType();
            cellTypes = type*ones(nelems,1);
            headCellTypes = '\nCELL_TYPES %d\n';
            fprintf(obj.outputFile, headCellTypes, nelems);
            fprintf(obj.outputFile, '%d ', cellTypes);
        end
        
        function writeData(obj)
            headPointData = '\n\nPOINT_DATA %d\n';
            headScalars = '\nSCALARS d_u float 3';
            headLookupT = '\nLOOKUP_TABLE default\n';
            nnodes = size(obj.coord,1);
%             ndofs = size(obj.d_u, 1);
            fprintf(obj.outputFile, headPointData, nnodes);
            fprintf(obj.outputFile, headScalars);
            fprintf(obj.outputFile, headLookupT);
            fprintf(obj.outputFile, '%f ', obj.d_u');
        end
        
        function saveFile(obj)
        end
        
        function t = getCellType(obj)
            switch obj.mesh.type
                case 'TRIANGLE'
                    t = 5;
                case 'QUAD'
                    t = 9;
                case 'TETRAHEDRA'
                    t = 10;
                case 'HEXAHEDRA'
                    t = 12;
            end
            
        end
        
    end
    
end

