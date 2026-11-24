function dataset = generateXiLatticeDataset(settings)
    if nargin < 1
        settings = latticeXiSettings();
    end
    ds = settings.dataset;
    if ~isfolder(ds.folder)
        mkdir(ds.folder);
    end
    generator = XiLatticeDatasetGenerator(settings);
    geometry  = buildGeometry(ds);
    scales    = generator.computeScales(ds.rho);
    checkManifest(settings,geometry,scales);
    pending   = findPending(ds.folder,size(geometry,1));
    pending   = sortByCost(generator,geometry,pending);
    fprintf('%d geometries in total, %d pending, %d workers\n',size(geometry,1),numel(pending),ds.nWorkers);
    runPending(settings,geometry,scales,pending);
    dataset = assembleDataset(settings,generator,geometry);
    save(fullfile(ds.folder,ds.outputFile),'-struct','dataset');
    fprintf('Dataset with %d samples saved in %s\n',size(dataset.data,1),fullfile(ds.folder,ds.outputFile));
end

function geometry = buildGeometry(ds)
    [A,T]    = ndgrid(ds.a(:),ds.t(:));
    geometry = [A(:),T(:)];
end

function checkManifest(settings,geometry,scales)
    file     = fullfile(settings.dataset.folder,'manifest.mat');
    manifest = struct('geometry',geometry,'scales',scales,'rho',settings.dataset.rho(:), ...
                      'material',settings.material,'mesh',settings.mesh, ...
                      'motif',settings.motif,'hole',settings.hole);
    if isfile(file)
        previous = load(file);
        if ~isequaln(previous.manifest,manifest)
            error('generateXiLatticeDataset:DifferentSettings', ...
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

function pending = sortByCost(generator,geometry,pending)
    cost = arrayfun(@(i) prod(generator.computeDivisions(geometry(i,:))),pending);
    [~,order] = sort(cost,'descend');
    pending   = pending(order);
end

function runPending(settings,geometry,scales,pending)
    nWorkers = settings.dataset.nWorkers;
    if nWorkers > 0 && isempty(gcp('nocreate'))
        parpool(nWorkers);
    end
    parfor (k = 1:numel(pending),nWorkers)
        iGeom = pending(k);
        computeAndSave(settings,geometry(iGeom,:),scales,iGeom);
    end
end

function computeAndSave(settings,params,scales,iGeom)
    tStart    = tic;
    generator = XiLatticeDatasetGenerator(settings);
    [C,rho,xi] = generator.computeHomogenizationAtScale(scales,params);
    result    = struct('params',params,'xi',xi,'C',C,'rho',rho,'rhoNominal',settings.dataset.rho(:),'time',toc(tStart));
    file      = geometryFile(settings.dataset.folder,iGeom);
    temporary = [file,'.tmp.mat'];
    save(temporary,'-struct','result');
    movefile(temporary,file,'f');
    fprintf('geometry %4d | a %.3f | t %.3f | %7.1f s\n',iGeom,params(1),params(2),result.time);
end

function file = geometryFile(folder,iGeom)
    file = fullfile(folder,sprintf('geometry_%05d.mat',iGeom));
end

function dataset = assembleDataset(settings,generator,geometry)
    symmetries = cellfun(@(name) generator.getSymmetry(name),settings.dataset.symmetries,'UniformOutput',false);
    blocks     = {};
    for iGeom = 1:size(geometry,1)
        result = load(geometryFile(settings.dataset.folder,iGeom));
        copies = expandBySymmetry(result.params,result.C,symmetries);
        for k = 1:numel(copies)
            blocks{end+1} = createBlock(generator,copies(k),result); %#ok<AGROW>
        end
    end
    rows   = concatenateBlocks([blocks{:}]);
    [~,order] = sortrows([rows.params,rows.rho]);
    dataset.data       = [rows.rho(order),rows.params(order,1),rows.xi(order),rows.values(order,:)];
    dataset.columns    = ["rho","a","xi",generator.getComponentNames()];
    dataset.t          = rows.params(order,2);
    dataset.tensors    = rows.tensors(:,:,:,:,order);
    dataset.rhoNominal = rows.rhoNominal(order);
    dataset.source     = rows.source(order);
    dataset.settings   = settings;
end

function rows = concatenateBlocks(blocks)
    rows.rho        = vertcat(blocks.rho);
    rows.rhoNominal = vertcat(blocks.rhoNominal);
    rows.params     = vertcat(blocks.params);
    rows.xi         = vertcat(blocks.xi);
    rows.values     = vertcat(blocks.values);
    rows.tensors    = cat(5,blocks.tensors);
    rows.source     = vertcat(blocks.source);
end

function copies = expandBySymmetry(params,C,symmetries)
    copies = struct('params',params,'C',C,'source',"computed");
    for iSym = 1:numel(symmetries)
        s = symmetries{iSym};
        for k = 1:numel(copies)
            if ~s.isInvariant(copies(k).params)
                copies(end+1) = struct('params',s.mapParams(copies(k).params), ...
                                       'C',s.mapTensor(copies(k).C), ...
                                       'source',copies(k).source + "+" + s.name); %#ok<AGROW>
            end
        end
    end
end

function block = createBlock(generator,copy,result)
    nRho             = numel(result.rho);
    block.rho        = result.rho(:);
    block.rhoNominal = result.rhoNominal(:);
    block.params     = repmat(copy.params,nRho,1);
    block.xi         = result.xi(:);
    block.values     = generator.tensorToVector(copy.C);
    block.tensors    = copy.C;
    block.source     = repmat(copy.source,nRho,1);
end