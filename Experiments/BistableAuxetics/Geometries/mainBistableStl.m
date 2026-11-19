here   = fileparts(mfilename('fullpath'));
outDir = fullfile(here,'..','results','stl');
if ~exist(outDir,'dir')
    mkdir(outDir);
end

specimens = {'hexagon',@PinwheelSpecimenGeometry.hexagonParams,2.3,4;
             'squareThick',@(n) PinwheelSpecimenGeometry.squareParams(n,0.05),2.3,2;
             'squareThin',@(n) PinwheelSpecimenGeometry.squareParams(n,0.025),2.3,2};
variants  = {'paper','print'};

for iS = 1:size(specimens,1)
    for nCells = [1 3]
        geo = PinwheelSpecimenGeometry(specimens{iS,2}(nCells));
        for iV = 1:2
            if iV == 1
                s.scale     = 1;
                s.thickness = specimens{iS,3};
            else
                s.scale     = specimens{iS,4};
                s.thickness = 3;
            end
            name = sprintf('%s_%dx%d_%s',specimens{iS,1},nCells,nCells,variants{iV});
            s.shape    = geo.shape;
            s.fileName = fullfile(outDir,[name,'.stl']);
            info = PinwheelStlExporter(s).export();
            reread = stlread(s.fileName);
            fprintf(['%s: size %.1f x %.1f x %.1f mm, faces %d (read back %d), closed %d, ', ...
                'volume %.1f, expected %.1f\n'],name,info.size,info.nFaces, ...
                size(reread.ConnectivityList,1),info.isClosed,info.volume,info.expected);
        end
    end
end
