classdef LatticeDatasetGenerator < handle

    properties (Constant, Access = private)
        outputComponents      = [1 1 1 1; 2 2 2 2; 1 1 2 2; 1 1 1 2; 2 2 1 2; 1 2 1 2]
        homogeneousDivisions  = 4
    end

    properties (Access = private)
        transformation
        meshSettings
        referenceCell
        homogeneousCell
        motif
        homogenizer
    end

    methods (Access = public)

        function obj = LatticeDatasetGenerator(cParams)
            obj.transformation = CellTransformation(cParams.transformation);
            obj.meshSettings   = cParams.mesh;
            obj.referenceCell  = obj.createCell(eye(2));
            obj.homogeneousCell = obj.createHomogeneousCell();
            obj.motif          = CellMotif(cParams.motif);
            obj.homogenizer    = MicroHomogenizer(cParams.material);
        end

        function C = computeHomogenizationAt(obj,rho,varargin)
            C = obj.computeHomogenizationAtScale(obj.computeScales(rho),varargin{:});
        end

        function [C,rho] = computeHomogenizationAtScale(obj,scale,varargin)
            T       = obj.transformation.compute([varargin{:}]);
            [C,rho] = obj.computeForMap(T,scale);
        end

        function [C,rho] = computeHomogenizationAtMap(obj,rho,T)
            [C,rho] = obj.computeForMap(T,obj.computeScales(rho));
        end

        function rho = computeVolumeFraction(obj,scale)
            mesh = obj.referenceCell.mesh;
            dens = obj.motif.computeDensity(mesh,eye(mesh.ndim),scale);
            rho  = obj.motif.computeVolumeFraction(mesh,dens);
        end

        function scales = computeScales(obj,rhoValues)
            referenceMesh = obj.referenceCell.mesh;
            scales = arrayfun(@(rho) obj.motif.computeScale(referenceMesh,rho),rhoValues(:));
        end

        function div = computeDivisions(obj,varargin)
            T       = obj.transformation.compute([varargin{:}]);
            [~,div] = PeriodicCell.discretize(T,obj.meshSettings.divUnit,obj.meshSettings.strategy);
        end

        function [X,Y,isSolid,T] = sampleMicrostructure(obj,scale,params,nCells,resolution)
            T      = obj.transformation.compute(params);
            xi     = linspace(0,nCells,resolution*nCells+1);
            [S1,S2] = ndgrid(xi,xi);
            x      = T*[S1(:)';S2(:)'];
            phi    = obj.motif.evaluateLevelSet(x,T,scale);
            X      = reshape(x(1,:),size(S1));
            Y      = reshape(x(2,:),size(S1));
            isSolid = reshape(phi > 0,size(S1));
        end

        function sym = getSymmetry(obj,name)
            sym = obj.transformation.getSymmetry(name);
        end

        function y = tensorToVector(obj,C)
            nSamples = size(C,5);
            nComp    = size(obj.outputComponents,1);
            y        = zeros(nSamples,nComp);
            for iComp = 1:nComp
                idx = num2cell(obj.outputComponents(iComp,:));
                y(:,iComp) = reshape(C(idx{:},:),[],1);
            end
        end

        function names = getComponentNames(obj)
            names = "C" + string(join(string(obj.outputComponents),"",2))';
        end

        function dataset = computeDataset(obj,rhoValues,geometryValues)
            scales   = obj.computeScales(rhoValues);
            nRho     = numel(scales);
            nGeom    = size(geometryValues,1);
            nColumns = 1 + size(geometryValues,2) + size(obj.outputComponents,1);
            dataset  = zeros(nRho*nGeom,nColumns);
            for iGeom = 1:nGeom
                params  = geometryValues(iGeom,:);
                [C,rho] = obj.computeForMap(obj.transformation.compute(params),scales);
                rows    = (iGeom-1)*nRho + (1:nRho);
                dataset(rows,:) = [rho,repmat(params,nRho,1),obj.tensorToVector(C)];
            end
        end

        function Chomog = computeGrid(obj,rhoValues,varargin)
            grids      = cell(1,numel(varargin));
            [grids{:}] = ndgrid(varargin{:});
            params     = cell2mat(cellfun(@(g) g(:),grids,'UniformOutput',false));
            scales     = obj.computeScales(rhoValues);
            tensorSize = repmat(obj.referenceCell.mesh.ndim,1,4);
            nGeom      = size(params,1);
            Chomog     = zeros([tensorSize,numel(scales),nGeom]);
            for iGeom = 1:nGeom
                T = obj.transformation.compute(params(iGeom,:));
                Chomog(:,:,:,:,:,iGeom) = obj.computeForMap(T,scales);
            end
            Chomog = reshape(Chomog,[tensorSize,numel(scales),cellfun(@numel,varargin)]);
        end

        function plotMicrostructures(obj,rhoValues,varargin)
            params       = [varargin{:}];
            T            = obj.transformation.compute(params);
            periodicCell = obj.createCell(T);
            scales       = obj.computeScales(rhoValues);
            for i = 1:numel(scales)
                dens = obj.motif.computeDensity(periodicCell.mesh,T,scales(i));
                figure
                dens.plot
                shading interp
                colormap(flipud(pink))
                axis equal
                title(sprintf('\\rho = %g, params = %s',rhoValues(i),mat2str(params,4)))
            end
        end

    end

    methods (Access = private)

        function c = createCell(obj,T)
            s   = obj.meshSettings;
            s.T = T;
            c   = PeriodicCell(s);
        end

        function c = createHomogeneousCell(obj)
            s.T        = eye(2);
            s.divUnit  = obj.homogeneousDivisions;
            s.strategy = 'Mapped';
            c = PeriodicCell(s);
        end

        function C = computeHomogeneous(obj,scale)
            mesh = obj.homogeneousCell.mesh;
            dens = obj.motif.computeDensity(mesh,eye(mesh.ndim),scale);
            C    = obj.homogenizer.compute(obj.homogeneousCell,dens);
        end

        function [C,rho] = computeForMap(obj,T,scales)
            periodicCell = obj.createCell(T);
            mesh         = periodicCell.mesh;
            nScales      = numel(scales);
            C            = zeros([repmat(mesh.ndim,1,4),nScales]);
            rho          = zeros(nScales,1);
            for i = 1:nScales
                if scales(i) == 0 || isinf(scales(i))
                    rho(i) = double(scales(i) == 0);
                    C(:,:,:,:,i) = obj.computeHomogeneous(scales(i));
                else
                    dens   = obj.motif.computeDensity(mesh,T,scales(i));
                    rho(i) = obj.motif.computeVolumeFraction(mesh,dens);
                    C(:,:,:,:,i) = obj.homogenizer.compute(periodicCell,dens);
                end
            end
        end

    end

end