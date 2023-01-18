classdef FunctionPrinter < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        
    end
    
    properties (Access = private)
        fun
        mesh
        filename
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
            d = obj.createMeshPostProcessDataBase('hey');
            p = MeshPrinter();
            p.print(d);
        end

        function d = createMeshPostProcessDataBase(obj,fileName)
            dI.mesh        = obj.mesh;
            dI.outFileName = fileName;
            ps = PostProcessDataBaseCreator(dI);
            db = ps.create();
            d = rmfield(db, 'gtype');
            d.ndim = size(obj.mesh.coord,2); % why not
        end

        function printResults(obj)
            fid = fopen('hey.flavia.res', 'w');
            obj.printResHeader(fid);
            obj.printResGaussInfo(fid);
            obj.printResFValues(fid);
            fclose(fid);
        end

        function printResHeader(obj, fid)
            fprintf(fid, 'GiD Post Results File 1.0 \n');
        end

        function printResGaussInfo(obj, fid)
            fprintf(fid, 'GaussPoints "Guass up?" Elemtype Triangle \n');
            fprintf(fid, 'Number of Gauss Points: 1 \n');
            fprintf(fid, 'Nodes not included \n');
            fprintf(fid, 'Natural Coordinates: given \n');
            fprintf(fid, ' 3.33333e-01  3.33333e-01  \n');
            fprintf(fid, 'End GaussPoints  \n');
        end

        function printResFValues(obj, fid)
            nodStr  = '%d ';
            numStr  = '%0.5g ';
            fValues = obj.fun.fValues;
            nNodes  = length(fValues);
            nodeMat = (1:nNodes)';
            results = [nodeMat, fValues];
            fValsStr  = repmat(numStr,[1, obj.fun.ndimf]);
            formatStr = [nodStr,fValsStr,'\n'];
            funTypeStr = obj.getFunctionTypeString();
            resHeaderStr = ['\nResult "obladi" "FunResults" 0', funTypeStr, 'OnNodes \n'];
            compsStr = obj.getComponentString();
            fprintf(fid, resHeaderStr);
            fprintf(fid, compsStr);
            fprintf(fid, 'Values  \n');
            fprintf(fid, formatStr,results');
            fprintf(fid, 'End Values\n\n');
            
        end

        function s = getFunctionTypeString(obj)
            if obj.fun.ndimf == 1
                s = ' Scalar ';
            else
                s = ' Vector ';
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