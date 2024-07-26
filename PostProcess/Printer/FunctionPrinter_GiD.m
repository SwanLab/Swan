classdef FunctionPrinter_GiD < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        
    end
    
    properties (Access = private)
        fun
        funNames
        mesh
        filename
        quadrature
    end
    
    methods (Access = public)
        
        function obj = FunctionPrinter_GiD(cParams)
            obj.init(cParams)
        end

        function appendFunction(obj, fun, name)
            obj.fun{end+1} = fun;
            obj.funNames{end+1} = name;
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
            obj.initFunctionNames(cParams);
        end
        
        function initFunctionNames(obj, cParams)
            if numel(obj.fun) == 1
                obj.funNames{1} = 'fValues';
            else
                obj.funNames = cParams.funNames;
            end
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
            for iFun = 1:numel(obj.fun)
                obj.printResFValues(fid, iFun);
            end
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

        function printResFValues(obj, fid, iFun)
            funTypeStr = obj.getFunctionTypeString(iFun);
            locTypeStr = obj.getFValuesLocationString(iFun);
            funNameStr = obj.funNames{iFun};
            resHeaderStr = ['\nResult "', funNameStr,'" "FunResults" 0', funTypeStr, locTypeStr,'\n'];
            compsStr = obj.getComponentString(iFun);
            [results,frmat] =  obj.fun{iFun}.getDataToPrint();
            fprintf(fid, resHeaderStr);
            fprintf(fid, compsStr);
            fprintf(fid, 'Values  \n');
            fprintf(fid,frmat,results{:});
            fprintf(fid, 'End Values\n\n');
            
        end

        function s = getFunctionTypeString(obj, iFun)
            if obj.fun{iFun}.ndimf == 1
                s = ' Scalar ';
            else
                s = ' Vector ';
            end
        end

        function s = getFValuesLocationString(obj, iFun)
            f  = obj.fun{iFun};
            type = class(f);
            switch type
                case 'LagrangianFunction'
                    switch f.order
                        case 'P0'
                            s = 'OnGaussPoints "Gauss Points" ';
                        otherwise
                            s = 'OnNodes ';
                    end
                case 'FGaussDiscontinuousFunction'
                    s = 'OnGaussPoints "Gauss Points" ';
            end
        end

        function s = getComponentString(obj, iFun)
            switch obj.fun{iFun}.ndimf
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