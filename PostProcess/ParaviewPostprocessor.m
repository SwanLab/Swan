classdef ParaviewPostprocessor < handle
    
    % Its output is a .vtu file. It can represent both simple fields (and
    % supposedly data at Gaussian points).
    
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
        
        function obj = ParaviewPostprocessor(cParams)
            obj.init(cParams);
            obj.openFile();
            obj.createPiece();
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
            fullfile = strcat(obj.filename, '.vtu');
            obj.outputFile = fopen(fullfile,'w');
        end
        
        function createPiece(obj)
            docNode = com.mathworks.xml.XMLUtils.createDocument('VTKFile');
            % Creating VTKFile
            fileN = docNode.getDocumentElement;
            fileN.setAttribute('type', 'UnstructuredGrid');
            fileN.setAttribute('version', '0.1');
            fileN.setAttribute('byte_order', 'LittleEndian');

            % Creating UnstructuredGrid
            gridN = docNode.createElement('UnstructuredGrid');
            fileN.appendChild(gridN);

            % Creating Piece
            pieceN = docNode.createElement('Piece');
            pieceN.setAttribute('NumberOfPoints', string(size(obj.coord,1))); %
            pieceN.setAttribute('NumberOfCells', string(size(obj.connec,1)));
            gridN.appendChild(pieceN);

            % Creating Points
            pointsN = obj.createPointsNode(docNode);
            pieceN.appendChild(pointsN);

            % Creating Cells
            cellsN = obj.createCellsNode(docNode);
            pieceN.appendChild(cellsN);

            % Create PointData
            pointDataN = docNode.createElement('PointData');
            pieceN.appendChild(pointDataN);
            
            % Create Displacement DataArray
            displacementDAN = obj.createDisplacementNode(docNode);
            pointDataN.appendChild(displacementDAN);
            
            text = xmlwrite(docNode);
            fprintf(obj.outputFile, text);
        end
        
        function n = createPointsNode(obj, docNode)
            n = docNode.createElement('Points');
            coordStr = obj.writeCoords();
            dan = docNode.createElement('DataArray');
            dan.setAttribute('type', 'Float32');
            dan.setAttribute('Name', 'Points');
            dan.setAttribute('NumberOfComponents', '3'); % !!
            dan.setAttribute('format', 'ascii');
            dat = docNode.createTextNode(coordStr);
            dan.appendChild(dat);
            n.appendChild(dan);
        end
        
        function n = createCellsNode(obj, docNode)
            n = docNode.createElement('Cells');
            [connecStr, offsetStr, typesStr] = obj.writeConnec();
            connecN = obj.createCellsConnecNode(docNode, connecStr);
            offset = obj.createCellsOffsetNode(docNode, offsetStr);
            types = obj.createCellsTypeNode(docNode, typesStr);
            n.appendChild(connecN);
            n.appendChild(offset);
            n.appendChild(types);
        end
        
        function n = createCellsConnecNode(obj, docNode, cstr)
            n = docNode.createElement('DataArray');
            n.setAttribute('type', 'Int32');
            n.setAttribute('Name', 'connectivity');
            n.setAttribute('format', 'ascii');
            t = docNode.createTextNode(cstr);
            n.appendChild(t);
        end
        
        function n = createCellsOffsetNode(obj, docNode, ostr)
            n = docNode.createElement('DataArray');
            n.setAttribute('type', 'Int32');
            n.setAttribute('Name', 'offsets');
            n.setAttribute('format', 'ascii');
            t = docNode.createTextNode(ostr);
            n.appendChild(t);
        end
        
        function n = createCellsTypeNode(obj, docNode, tstr)
            n = docNode.createElement('DataArray');
            n.setAttribute('type', 'UInt8');
            n.setAttribute('Name', 'types');
            n.setAttribute('format', 'ascii');
            t = docNode.createTextNode(tstr);
            n.appendChild(t);
        end

        function n = createDisplacementNode(obj, docNode)
            dispStr = sprintf('%.4f %.4f %.4f\n', obj.d_u');
            n = docNode.createElement('DataArray');
            n.setAttribute('type', 'Float64');
            n.setAttribute('Name', 'Displacement');
            n.setAttribute('NumberOfComponents', '3'); % !!
            n.setAttribute('format', 'ascii');
            t = docNode.createTextNode(dispStr);
            n.appendChild(t);
        end
        
        function coordStr = writeCoords(obj)
            nnodes = size(obj.coord,1);
            if (size(obj.coord,2) == 2)
                obj.coord = [obj.coord, zeros(nnodes,1)];
            end
            coordStr = sprintf('%.4f %.4f %.4f\n', obj.coord');
        end
        
        function [connecStr, offsetStr, typesStr] = writeConnec(obj)
            nelems = size(obj.connec,1);
            nnodEl = size(obj.connec,2);
            connecP = obj.connec-1;
            format = [repmat(' %d',1,nnodEl),'\n'];
            connecStr = sprintf(format, connecP');
            
            offset = nnodEl:nnodEl:nelems*nnodEl;
            offsetStr = sprintf('%d ', offset);
            
            vtkType = obj.getCellType();
            types = vtkType*ones(nelems,1);
            typesStr = sprintf('%d ', types);
        end
        
        function writeData(obj)
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

