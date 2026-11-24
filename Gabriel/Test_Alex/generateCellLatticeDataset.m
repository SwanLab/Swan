function dataset = generateCellLatticeDataset(settings)
    ds = settings.dataset;
    if ~isfolder(ds.folder)
        mkdir(ds.folder);
    end
    generator = SidesLatticeDatasetGenerator(settings);
    geometry  = buildGeometry(ds);
    checkManifest(settings,geometry);
    pending   = findPending(ds.folder,size(geometry,1));
    pending   = sortByCost(generator,ds,geometry,pending);
    fprintf('%d geometries in total, %d pending, %d workers\n',size(geometry,1),numel(pending),ds.nWorkers);
    runPending(settings,geometry,pending);
    dataset = assembleDataset(settings,generator,geometry);
    save(fullfile(ds.folder,ds.outputFile),'-struct','dataset');
    fprintf('Dataset with %d samples saved in %s\n',size(dataset.data,1),fullfile(ds.folder,ds.outputFile));
end

function geometry = buildGeometry(ds)
    [iCell,M1] = ndgrid(1:size(ds.cells,1),ds.m(:));
    geometry   = [iCell(:),M1(:)];
end

function sides = buildSides(ds,m1)
    sides = [m1*ones(numel(ds.m),1),ds.m(:)];
end

function checkManifest(settings,geometry)
    file     = fullfile(settings.dataset.folder,'manifest.mat');
    manifest = struct('geometry',geometry,'cells',settings.dataset.cells,'m',settings.dataset.m(:), ...
                      'material',settings.material,'mesh',settings.mesh, ...
                      'motif',settings.motif,'hole',settings.hole);
    if isfile(file)
        previous = load(file);
        if ~isequaln(previous.manifest,manifest)
            error('generateCellLatticeDataset:DifferentSettings', ...
                'The folder %s contains a dataset generated with different settings. Use another folder.',settings.dataset.folder);
        end
    else
        save(file,'manifest');
    end
end

function pending = findPending(folder,nGeom)
    done    = arrayfun(@(i) isfile(geometryFile(folder,i)),1:nGeom);
    pending = find(~done);
end

function pending = sortByCost(generator,ds,geometry,pending)
    cost = arrayfun(@(i) prod(generator.computeDivisions(ds.cells(geometry(i,1),:))),pending);
    [~,order] = sort(cost,'descend');
    pending   = pending(order);
end

function runPending(settings,geometry,pending)
    nWorkers = settings.dataset.nWorkers;
    if nWorkers > 0 && isempty(gcp('nocreate'))
        parpool(nWorkers);
    end
    parfor (k = 1:numel(pending),nWorkers)
        iGeom = pending(k);
        computeAndSave(settings,geometry(iGeom,:),iGeom);
    end
end

function computeAndSave(settings,params,iGeom)
    tStart    = tic;
    generator = SidesLatticeDatasetGenerator(settings);
    cellParams = settings.dataset.cells(params(1),:);
    sides     = buildSides(settings.dataset,params(2));
    [C,rho]   = generator.computeHomogenizationAt(cellParams,sides);
    result    = struct('cellParams',cellParams,'sides',sides,'C',C,'rho',rho,'time',toc(tStart));
    file      = geometryFile(settings.dataset.folder,iGeom);
    temporary = [file,'.tmp.mat'];
    save(temporary,'-struct','result');
    movefile(temporary,file,'f');
    fprintf('geometry %5d | a %7.3f | s %5.2f | r %5.3f | m1 %.3f | %7.1f s\n',iGeom,cellParams(1),cellParams(2),cellParams(3),params(2),result.time);
end

function file = geometryFile(folder,iGeom)
    file = fullfile(folder,sprintf('geometry_%05d.mat',iGeom));
end

function dataset = assembleDataset(settings,generator,geometry)
    ds             = settings.dataset;
    transformation = createLatticeTransformation(settings);
    computedShapes = unique(round(ds.cells(:,2:3),10),'rows');
    nGeom          = size(geometry,1);
    params         = cell(nGeom,2);
    rho            = cell(nGeom,2);
    tensors        = cell(nGeom,2);
    source         = cell(nGeom,2);
    for iGeom = 1:nGeom
        result = load(geometryFile(ds.folder,iGeom));
        n      = size(result.sides,1);
        params{iGeom,1}  = [result.sides,repmat(result.cellParams,n,1)];
        rho{iGeom,1}     = result.rho(:);
        tensors{iGeom,1} = result.C;
        source{iGeom,1}  = repmat("computed",n,1);
        twin = transformation.computeTwin(params{iGeom,1});
        if ~ismember(round(twin(1,4:5),10),computedShapes,'rows')
            params{iGeom,2}  = twin;
            rho{iGeom,2}     = result.rho(:);
            tensors{iGeom,2} = result.C;
            source{iGeom,2}  = repmat("twin",n,1);
        end
    end
    allParams  = vertcat(params{:});
    allTensors = cat(5,tensors{:});
    allRho     = vertcat(rho{:});
    allSource  = vertcat(source{:});
    [~,order]  = sortrows(round(allParams(:,[4,5,3,1,2]),10));
    allParams  = allParams(order,:);
    allTensors = allTensors(:,:,:,:,order);
    dataset.data       = [allParams,allRho(order),generator.tensorToVector(allTensors)];
    dataset.columns    = ["m1","m2","a","s","r","rho",generator.getComponentNames()];
    dataset.rhoNominal = SidesLatticeDatasetGenerator.computeRho(allParams(:,1),allParams(:,2));
    dataset.xi         = SidesLatticeDatasetGenerator.computeXi(allParams(:,1),allParams(:,2));
    dataset.tensors    = allTensors;
    dataset.source     = allSource(order);
    dataset.settings   = settings;
end