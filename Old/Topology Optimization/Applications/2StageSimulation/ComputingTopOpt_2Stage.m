function ComputingTopOpt_2Stage
% checkout: git checkout origin/master -- path/to/file
% PENDING: AUTOMATIZE THIS FUNCTION WITH GID PICTUREsss
rho0Name = 'TFM.mat';
jumpTo2ndPart = false;

if jumpTo2ndPart == false

%     fileName = 'jaCantilever';
%     % Data input
%     s.testName = [fileName,'.m'];
%     s.x1       = 2;
%     s.y1       = 1;
%     s.N        = 150;
%     s.M        = 75;
%     s.P        = -100;
%     s.DoF      = 2;
% 
%     FEMWriter = FEMInputWriter(s);
%     FEMWriter.createTest;




    fileName = 'test_anisotropy_cantilever';
    s = Settings(fileName);
    s.warningHoleBC = false;
    s.printIncrementalIter = false;
    s.printChangingFilter = false;
    s.printing = false;
    translator = SettingsTranslator();
    translator.translate(s);
    fileName = translator.fileName;
    settings  = SettingsTopOptProblem(fileName);
    topOptSolver = TopOpt_Problem(settings);
    while topOptSolver.incrementalScheme.hasNext()
        topOptSolver.incrementalScheme.next();
        topOptSolver.optimizer.solveProblem();
    end

    rho0 = topOptSolver.designVariable.value;
    save(rho0Name,'rho0');

    % GiD Image Capturer:
%     outPutImageName = '"testgidpic"';
%     gidPath = '/home/joseantonio/GiDx64/gid-15.0.4/';
%     tclFile = 'callGiDCapturer.tcl';
%     tclFileTocall = 'CaptureImage3.tcl';
%     fid = fopen('/home/joseantonio/Documentos/GitHub/Swan/PostProcess/ImageCapturer/callGiDCapturer.tcl','w+');
%     fprintf(fid,['set path "/home/joseantonio/Documentos/GitHub/Swan/PostProcess/ImageCapturer','"\n']);
%     fprintf(fid,['set tclFile "',tclFileTocall,'"\n']);
%     fprintf(fid,['source $path$tclFile \n']);
%     fprintf(fid,['set output ',outPutImageName,' \n']);
%     fprintf(fid,['set inputFile "',fileName,'"\n']);
%     fprintf(fid,['CaptureImage $inputFile $output \n']);
%     fclose(fid);
%     command = [gidPath,'gid_offscreen -offscreen -t "source ',tclFile,'"'];
%     system(command);
%     inputImage  = [' ',outPutImageName,'.png'];
%     outPutImage = inputImage;
%     convert     = 'convert -crop 442x442+0+0 -gravity Center';
%     command = strcat(convert,' ',inputImage,' ',outPutImage);
%     system(command);
else
    load(rho0Name);
    fileName = 'test_anisotropy_cantilever_rho0';

    s = Settings(fileName);
    s.warningHoleBC = false;
    s.printIncrementalIter = false;
    s.printChangingFilter = false;
    s.printing = false;
    translator = SettingsTranslator();
    translator.translate(s);
    fileName = translator.fileName;
    settings  = SettingsTopOptProblem(fileName);
    DesignVariable = convertCharsToStrings(settings.designVarSettings.type);
    if  DesignVariable == "Density"
        settings.designVarSettings.creatorSettings.type = 'Given';
        settings.designVarSettings.creatorSettings.rho0 = rho0;
    elseif DesignVariable == "LevelSet"
        settings.designVarSettings.initialCase = 'given';
        settings.designVarSettings.creatorSettings.value = rho0;
    end
    topOptSolver = TopOpt_Problem(settings);
    while topOptSolver.incrementalScheme.hasNext()
        topOptSolver.incrementalScheme.next();
        topOptSolver.optimizer.solveProblem();
    end
end

end