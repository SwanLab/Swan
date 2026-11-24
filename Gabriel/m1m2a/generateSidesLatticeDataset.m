function dataset = generateSidesLatticeDataset(settings)
    if nargin < 1
        settings = latticeSidesSettings();
    end
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
    [A,M1]   = ndgrid(ds.a(:),ds.m(:));
    geometry = [A(:),M1(:)];
end

function sides = buildSides(ds,m1)
    % m2    = ds.m(ds.m <= m1);
    m2    = ds.m;
    if any(strcmp(ds.symmetries,'swapSides'))
        m2 = ds.m(ds.m <= m1);
    end
    sides = [m1*ones(numel(m2),1),m2(:)];
end

function checkManifest(settings,geometry)
    file     = fullfile(settings.dataset.folder,'manifest.mat');
    manifest = struct('geometry',geometry,'m',settings.dataset.m(:), ...
                      'material',settings.material,'mesh',settings.mesh, ...
                      'motif',settings.motif,'hole',settings.hole);
    if isfile(file)
        previous = load(file);
        if ~isequaln(previous.manifest,manifest)
            error('generateSidesLatticeDataset:DifferentSettings', ...
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
    cost = arrayfun(@(i) prod(generator.computeDivisions(geometry(i,1)))*size(buildSides(ds,geometry(i,2)),1),pending);
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
    a         = params(1);
    sides     = buildSides(settings.dataset,params(2));
    [C,rho]   = generator.computeHomogenizationAt(a,sides);
    result    = struct('a',a,'sides',sides,'C',C,'rho',rho,'time',toc(tStart));
    file      = geometryFile(settings.dataset.folder,iGeom);
    temporary = [file,'.tmp.mat'];
    save(temporary,'-struct','result');
    movefile(temporary,file,'f');
    fprintf('geometry %4d | a %.3f | m1 %.3f | %3d cells | %7.1f s\n',iGeom,a,params(2),size(sides,1),result.time);
end

function file = geometryFile(folder,iGeom)
    file = fullfile(folder,sprintf('geometry_%05d.mat',iGeom));
end

function dataset = assembleDataset(settings,generator,geometry)
    symmetries = cellfun(@(name) generator.getSymmetry(name),settings.dataset.symmetries,'UniformOutput',false);
    blocks     = {};
    for iGeom = 1:size(geometry,1)
        result = load(geometryFile(settings.dataset.folder,iGeom));
        params = [result.sides,result.a*ones(size(result.sides,1),1)];
        copies = expandBySymmetry(params,result.C,result.rho,symmetries);
        for k = 1:numel(copies)
            blocks{end+1} = createBlock(generator,copies(k)); %#ok<AGROW>
        end
    end
    rows      = concatenateBlocks([blocks{:}]);
    [~,order] = sortrows(rows.params);
    params    = rows.params(order,:);
    dataset.data       = [params,rows.rho(order),rows.values(order,:)];
    dataset.columns    = ["m1","m2","a","rho",generator.getComponentNames()];
    dataset.rhoNominal = SidesLatticeDatasetGenerator.computeRho(params(:,1),params(:,2));
    dataset.xi         = SidesLatticeDatasetGenerator.computeXi(params(:,1),params(:,2));
    dataset.tensors    = rows.tensors(:,:,:,:,order);
    dataset.source     = rows.source(order);
    dataset.settings   = settings;
end

function rows = concatenateBlocks(blocks)
    rows.params  = vertcat(blocks.params);
    rows.rho     = vertcat(blocks.rho);
    rows.values  = vertcat(blocks.values);
    rows.tensors = cat(5,blocks.tensors);
    rows.source  = vertcat(blocks.source);
end

function copies = expandBySymmetry(params,C,rho,symmetries)
    copies = struct('params',params,'C',C,'rho',rho,'source',repmat("computed",size(params,1),1));
    for iSym = 1:numel(symmetries)
        s = symmetries{iSym};
        for k = 1:numel(copies)
            keep = ~s.isInvariant(copies(k).params);
            if any(keep)
                copies(end+1) = struct('params',s.mapParams(copies(k).params(keep,:)), ...
                                       'C',s.mapTensor(copies(k).C(:,:,:,:,keep)), ...
                                       'rho',copies(k).rho(keep), ...
                                       'source',copies(k).source(keep) + "+" + s.name); %#ok<AGROW>
            end
        end
    end
end

function block = createBlock(generator,copy)
    block.params  = copy.params;
    block.rho     = copy.rho(:);
    block.values  = generator.tensorToVector(copy.C);
    block.tensors = copy.C;
    block.source  = copy.source(:);
end