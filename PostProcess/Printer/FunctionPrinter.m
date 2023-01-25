classdef FunctionPrinter < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        
    end
    
    properties (Access = private)
        fun
        mesh
        filename
        quadrature
    end
    
    methods (Access = public)
        
        function obj = FunctionPrinter(cParams)
            obj.init(cParams)
            
        end

        function print(obj)
            obj.printMesh();
            obj.printResults();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.fun      = cParams.fun;
            obj.mesh     = cParams.mesh;
            obj.filename = cParams.filename;
        end
        
        function printMesh(obj)
            d = obj.createMeshPostProcessDataBase(obj.filename);
            p = MeshPrinter();
            p.print(d);
        end

        function d = createMeshPostProcessDataBase(obj,fileName)
            dI.mesh        = obj.mesh;
            dI.outFileName = fileName;
            ps = PostProcessDataBaseCreator(dI);
            db = ps.create();
            d = rmfield(db, 'gtype');
            d.ndim = size(obj.mesh.coord,2); % yolo
        end

        function printResults(obj)
            f = obj.filename;
            fid = fopen([f, '.flavia.res'], 'w');
            obj.printResHeader(fid);
            obj.printResGaussInfo(fid);
            obj.printResFValues(fid);
            fclose(fid);
        end

        function printResHeader(obj, fid)
            fprintf(fid, 'GiD Post Results File 1.0 \n');
        end

        function printResGaussInfo(obj, fid)
            obj.createQuadrature();
            nDim  = size(obj.quadrature.posgp,1);
            nGaus = obj.quadrature.ngaus;
            el = obj.getGiDElementType();
            fprintf(fid, ['GaussPoints "Gauss Points" Elemtype ', el, '\n']);
            fprintf(fid, ['Number of Gauss Points: ', num2str(nGaus),' \n']);
            fprintf(fid, 'Nodes not included \n');
            fprintf(fid, 'Natural Coordinates: given \n');
            printFormat = [repmat('%12.5d ',1,nDim),'\n'];
            toPrint = obj.quadrature.posgp;
            fprintf(fid,printFormat,toPrint);
            fprintf(fid, 'End GaussPoints  \n');
        end
        
        function createQuadrature(obj)
            q = Quadrature.set(obj.mesh.type);
            q.computeQuadrature('LINEAR');
            obj.quadrature = q;
        end

        
        function s = getGiDElementType(obj)
            switch  obj.mesh.type
                case 'TRIANGLE'
                    s = 'Triangle';
                case 'QUAD'
                    s = 'Quadrilateral';
                case 'TETRAHEDRA'
                    s = 'Tetrahedra';
                case 'HEXAHEDRA'
                    s = 'Hexahedra';
            end
        end

        function printResFValues(obj, fid)
            nodStr  = '%d ';
            numStr  = '%0.5g ';
            [results,frmat] = obj.formatResultsMat();
            fValsStr  = repmat(numStr,[1, obj.fun.ndimf]);
            formatStr = [nodStr,fValsStr,'\n'];
            funTypeStr = obj.getFunctionTypeString();
            locTypeStr = obj.getFValuesLocationString();
            resHeaderStr = ['\nResult "fValues" "FunResults" 0', funTypeStr, locTypeStr,'\n'];
            compsStr = obj.getComponentString();
            fprintf(fid, resHeaderStr);
            fprintf(fid, compsStr);
            fprintf(fid, 'Values  \n');

            fprintf(fid,frmat,results{:});
%             fprintf(fid, formatStr,results');


%             pformat = ['%s ',repmat('%12.5d ',1,obj.fun.ndimf),'\n'];
%             c  = obj.computeElementStringColum();
%             fM = obj.computeTensorValueColums();
%             str = [c,fM]';
%             fprintf(fid,pformat,str{:});

            fprintf(fid, 'End Values\n\n');
            
        end

        
        function c = computeElementStringColum(obj)
            nElem = size(obj.mesh.connec, 1);
            nGaus = obj.quadrature.ngaus;
            allElem(:,1) = 1:nElem;
            colWidth = size(num2str(nElem),2);
            strInCol = repmat(' ',nElem*nGaus,colWidth);
            numIndex = 1:nGaus:nElem*nGaus;
            strInCol(numIndex,:) = num2str(allElem);
            c = cellstr(strInCol);
        end
        
        function fM = computeTensorValueColums(obj)
            fV = obj.fun.fValues;
            nGaus   = obj.quadrature.ngaus;
            nComp   = size(obj.fun.fValues, 2);
            nElem   = size(obj.mesh.connec, 1);
            fM  = zeros(nGaus*nElem,nComp);
            for istre = 1:nComp
                for igaus = 1:nGaus
                    rows = linspace(igaus,(nElem - 1)*nGaus + igaus,nElem);
                    fM(rows,istre) = fV(igaus,istre,:);
                end
            end
            fM = num2cell(fM);
        end

        function [res, format] = formatResultsMat(obj)
            % ngaus!! --> see GaussFieldPrinter
            % move to fun.getDataToPrint()
            switch class(obj.fun)
                case 'P1Function'
                    fValues = squeeze(obj.fun.fValues);
                    nNodes  = length(fValues);
                    nodeMat = (1:nNodes)';
                    res = [nodeMat, fValues];
                case 'P1DiscontinuousFunction'
                    ndims   = size(obj.fun.fValues, 1);
                    nelem   = size(obj.mesh.connec, 1);
                    nnodeEl = size(obj.mesh.connec, 2);
                    fV = reshape(obj.fun.fValues, [ndims, nelem*nnodeEl])';
                    nNodes  = length(fV);
                    nodeMat = (1:nNodes)';
                    res = [nodeMat, fV];
                case 'P0Function'
                    [res, format] = obj.fun.getDataToPrint();
                case 'FGaussDiscontinuousFunction'
                    [res, format] = obj.fun.getDataToPrint();

            end
        end

        function s = getFunctionTypeString(obj)
            if obj.fun.ndimf == 1
                s = ' Scalar ';
            else
                s = ' Vector ';
            end
        end

        function s = getFValuesLocationString(obj)
            switch class(obj.fun)
                case {'P1Function', 'P1DiscontinuousFunction'}
                    s = 'OnNodes ';
                case {'P0Function', 'FGaussDiscontinuousFunction'}
                    s = 'OnGaussPoints "Gauss Points" ';
            end
        end

        function s = getComponentString(obj)
            switch obj.fun.ndimf
                case 1
                    s = 'ComponentNames  "x"\n';
                case 2
                    s = 'ComponentNames  "x1", "x2"\n';
                case 3
                    s = 'ComponentNames  "x1", "x2", "x3"\n';
                case 4
                    s = 'ComponentNames  "x1", "x2", "x3", "x4"\n';
                case 6
                    s = 'ComponentNames  "x1", "x2", "x3", "x4", "x5", "x6"\n';
            end
        end

    end
    
end