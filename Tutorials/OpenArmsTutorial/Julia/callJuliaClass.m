function output = callJuliaClass(className, methodName, params)
    % Path to Julia scripts folder  
    currentScriptPath = fileparts(mfilename('fullpath'));
    juliaScriptDir = fullfile(currentScriptPath, className);

    % Create filenames
    inputFile = 'julia_input.json';
    outputFile = 'julia_output.json';

    % Add output field so Julia knows where to write
    params.output = outputFile;

    % Write input JSON
    fid = fopen(inputFile, 'w');
    fwrite(fid, jsonencode(params));
    fclose(fid);

    % Build Julia command
    juliaScriptName = ['call_' className '_' methodName '.jl'];
    juliaScript = fullfile(juliaScriptDir, juliaScriptName);
    command = ['julia ' juliaScript ' ' inputFile];
    system(command);

    % Read output JSON
    fid = fopen(outputFile, 'r');
    raw = fread(fid, '*char')';
    fclose(fid);

    output = jsondecode(raw);
end