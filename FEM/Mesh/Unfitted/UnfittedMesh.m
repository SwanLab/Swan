classdef UnfittedMesh < handle
    
    properties (GetAccess = public, SetAccess = private)
        innerMesh
        innerCutMesh
        boundaryCutMesh
        unfittedBoundaryMesh
        backgroundMesh
    end
    
    properties (Access = private)
        cutMesh
        fullCells
        cutCells
        emptyCells
        plotter
    end
    
    properties (Access = private)
        boundaryMesh
        isInBoundary
        levelSet
    end
    
    methods (Access = public)
        
        function obj = UnfittedMesh(cParams)
            obj.init(cParams);
        end
        
        function compute(obj,lSet)
            obj.levelSet = lSet;
            obj.classifyCells();
            obj.computeInnerMesh();
            if ~isempty(obj.cutCells)
                obj.computeCutMesh();
                obj.computeInnerCutMesh();
                obj.computeBoundaryCutMesh();
            else 
                obj.cutMesh         = [];
                obj.innerCutMesh    = [];
                obj.boundaryCutMesh = [];
            end
            obj.computeUnfittedBoxMesh();
            obj.createPlotter();
        end
        
        function createPlotter(obj)
            s.uMesh = obj;
            obj.plotter = UnfittedMeshPlotter(s);
        end
        
        function plotBoundary(obj)
            obj.plotter.plotBoundary();
        end
        
        function plotStructureInColor(obj,color)
            obj.plotter.plotDomainInColor(color);
        end
        
        function plot(obj)
            obj.plotter.plotDomain();
        end
        
        function plotAll(obj)
            obj.plotter.plotAll();
        end
        
        function plotNormals(obj)
            if ~isempty(obj.boundaryCutMesh)
                obj.boundaryCutMesh.mesh.plotNormals();
            end
        end
        
        function dv = computeDvolume(obj,quad)
            if ~isempty(obj.innerMesh)
                dvI = obj.computeInnerDvolume(quad);
            else 
                dvI = 0;
            end
            dvC = obj.computeInnerCutDvolume(quad);
            dv = dvC + dvI;
        end
        
        function print(obj, filename)
            d = obj.createPostProcessDataBase(filename);
            d.fields           = {obj.levelSet};
            d.printMode        = 'DesignVariable';
            d.nDesignVariables = 1;
            postProcess = Postprocess('TopOptProblem',d);
            postProcess.print(1,d)
        end
        
        function exportGiD(obj, cParams)
            % SAMPLE CPARAMS
            % s.filename        = 'hellothere';
            % s.resFilePath     = '/home/ton/test_micro23.flavia.res';
            % s.gidProjectPath  = '/home/ton/test_micro_project.gid';
            % s.meshElementSize = '0.0707107';
            % s.meshFileName    = 'hmmmm22';
            % s.swanPath        = '/home/ton/Github/Swan/';
            % s.gidPath         = '/home/ton/GiDx64/gid-16.1.2d/';
            filename = cParams.filename;
            obj.print(filename);
            cParams.resFilePath = [cParams.swanPath, 'Output/',filename,'/',filename,'1.flavia.res'];
            obj.createSurfaceMeshTclFromTemplate(cParams);
            pathTcl  = [cParams.swanPath,'PostProcess/STL/'];
            
            tclFile = [pathTcl,'CreateSurfaceMeshFile.tcl" '];
            command = [cParams.gidPath,'gid_offscreen -offscreen -t "source ',tclFile];
            system(command);
        end

        function m = exportInnerMesh(obj)
            s.unfittedMesh = obj;
            s.type = 'Matlab';
            s.filename= 'none'; % only for gid
            ime = InnerMeshExporter(s);
            m = ime.export();
        end

    end
    
    methods (Access = private)
        
        function dvolume = computeInnerDvolume(obj,quad)
            nelem = obj.backgroundMesh.nelem;
            ngaus = quad.ngaus;
            dv    = obj.innerMesh.mesh.computeDvolume(quad);
            iCells = obj.innerMesh.fullCells;
            dvolume = zeros(nelem,ngaus);
            dvolume(iCells,:) = dv;
        end
        
        function dvolume = computeInnerCutDvolume(obj,quad)
            nelem = obj.backgroundMesh.nelem;
            ngaus = quad.ngaus;
            dvolume = zeros(nelem,ngaus);
            if ~isempty(obj.innerCutMesh)
            cCells = obj.innerCutMesh.cellContainingSubcell;
            dv = obj.innerCutMesh.mesh.computeDvolume(quad);
            for igaus = 1:ngaus
               dvolume(:,igaus) = accumarray(cCells,dv(igaus,:)',[nelem,1],@sum,0);
            end
            end
        end
        
        function init(obj,cParams)
            obj.backgroundMesh = cParams.backgroundMesh;
            obj.boundaryMesh   = cParams.boundaryMesh;
        end
        
        function classifyCells(obj)
            allCells(:,1) = 1:obj.backgroundMesh.nelem;
            lsNodes  = obj.computelsNodes();
            isLsNeg  = lsNodes < 0;
            full  = all(isLsNeg,2);
            empty = all(~isLsNeg,2);
            cut   = ~or(full,empty);
            obj.fullCells  = allCells(full);
            obj.emptyCells = allCells(empty);
            obj.cutCells   = allCells(cut);
        end
        
        function ls = computelsNodes(obj)
            nodes = obj.backgroundMesh.connec;
            ls = zeros(size(nodes));
            for iNode = 1:size(nodes,2)
                ls(:,iNode) = obj.levelSet(nodes(:,iNode));
            end
        end
        
        function computeInnerMesh(obj)
            if ~isempty(obj.fullCells)
                s.backgroundMesh = obj.backgroundMesh;
                s.fullCells      = obj.fullCells;
                obj.innerMesh = InnerMesh(s);
            else
                obj.innerMesh = [];
            end
        end
        
        function computeCutMesh(obj)
            s.backgroundMesh = obj.backgroundMesh;
            s.cutCells       = obj.cutCells;
            s.levelSet       = obj.levelSet;
            obj.cutMesh = CutMesh.create(s);
            obj.cutMesh.compute();
        end
        
        function computeInnerCutMesh(obj)
            obj.innerCutMesh = obj.cutMesh.innerCutMesh;
        end
        
        function computeBoundaryCutMesh(obj)
            if ~isequal(obj.backgroundMesh.geometryType,'Line')
                obj.boundaryCutMesh = obj.cutMesh.boundaryCutMesh;
            end
            
        end
        
        function computeUnfittedBoxMesh(obj)
            s.boundaryMesh = obj.boundaryMesh;
            obj.unfittedBoundaryMesh = UnfittedBoundaryMesh(s);
            if ~isempty(obj.unfittedBoundaryMesh.meshes)
                ls   = obj.levelSet;
                obj.unfittedBoundaryMesh.compute(ls);
            end
        end

        function d = createPostProcessDataBase(obj, filename)
            d.mesh    = obj.backgroundMesh;
            d.outFileName = filename;
            d.ptype   = 'TopOpt';
            ps = PostProcessDataBaseCreator(d);
            d  = ps.create();
            d.ndim       = obj.backgroundMesh.ndim;
            d.pdim       = obj.backgroundMesh.ndim;
            d.designVar  = 'LevelSet';
        end

        function createSurfaceMeshTclFromTemplate(obj, cParams)
            resFilePath     = cParams.resFilePath;
            gidProjectPath  = cParams.gidProjectPath;
            meshElementSize = cParams.meshElementSize;
            meshFileName    = cParams.meshFileName;

            templateText = fileread('CreateSurfaceMeshFile_Template.tcl');
            targetFile = ['PostProcess/STL/', 'CreateSurfaceMeshFile.tcl'];
            fid = fopen(targetFile, 'w');
            fprintf(fid, ['set input_post_res "', resFilePath, '"\n']);
            fprintf(fid, ['set output_gid_project_name "', gidProjectPath, '"\n']);
            fprintf(fid, ['set mesh_element_size "', meshElementSize, '"\n']);
            fprintf(fid, ['set mesh_name "', meshFileName, '"\n']);
            fprintf(fid, ['\n',templateText]);
            fclose(fid);
        end
        
    end
    
    methods (Access = public)
        
        function mass = computeMass(obj)
            npnod = obj.backgroundMesh.nnodes;
            f = ones(npnod,1);
            s.mesh = obj;
            s.type = 'Unfitted';
            integrator = RHSintegrator.create(s);
            fInt = integrator.integrateInDomain(f);
            %%Now to check IntegrateNodal, later by obj.mesh.computeMass
            %disp('Interior')
            %sum(fInt<0)/size(fInt,1)
            mass = sum(fInt);
        end
        
        function mass = computePerimeter(obj)
            npnod = obj.backgroundMesh.nnodes;
            f = ones(npnod,1);
            s.mesh = obj;
            s.type = 'Unfitted';
            integrator = RHSintegrator.create(s);
            fInt = integrator.integrateInBoundary(f);
            %%Now to check IntegrateNodal, later by obj.mesh.computeMass
            %disp('Boundary')
            %sum(fInt<0)/size(fInt,1)
            mass = sum(fInt);
        end
        
    end
    
end