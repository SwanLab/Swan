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

        function printResults(obj)
            d = obj.createResultsPostProcessDataBase('hey');
            p = ResultsPrinter.create('FeFunction', d);
            p.print('',d);
        end
        
        function d = createResultsPostProcessDataBase(obj,fileName)
            quad = Quadrature.set(obj.mesh.type);
            quad.computeQuadrature('LINEAR');

            dI.mesh        = obj.mesh;
            dI.outFileName = fileName;
            ps = PostProcessDataBaseCreator(dI);
            d = ps.create();
            d.ndim = size(obj.mesh.coord,2);
            d.resultsDir = ''; % ep!!
            d.name = 'whatever';
            d.fields{1} = obj.fun.fValues;
            d.quad = quad; % !!
        end

        function d = createMeshPostProcessDataBase(obj,fileName)
            dI.mesh        = obj.mesh;
            dI.outFileName = fileName;
            ps = PostProcessDataBaseCreator(dI);
            db = ps.create();
            d = rmfield(db, 'gtype');
            d.ndim = size(obj.mesh.coord,2); % why not
        end

    end
    
end