clc; clear; close all;

% Set up path
addpath(genpath(pwd));
try
    rmpath(genpath('Old'));
catch
end

fid = fopen('TutorialShells.m', 'r');
code = fread(fid, '*char')';
fclose(fid);

code = strrep(code, 'properties (Access=private)', 'properties (Access=public)');
code = strrep(code, 'methods (Access = private)', 'methods (Access = public)');
% Ensure plotMatlab is 0 so it doesn't try to plot
code = strrep(code, 'plotMatlab = 1', 'plotMatlab = 0');
code = strrep(code, 'printParaview = 1', 'printParaview = 0');
code = strrep(code, 'problemType    = ''FORCED_VIBRATIONS'';', 'problemType = ''STATIC'';');

code50 = strrep(code, 'classdef TutorialShells < handle', 'classdef TutorialShells50 < handle');
code50 = strrep(code50, 'function obj = TutorialShells()', 'function obj = TutorialShells50()');
code50 = strrep(code50, 'elements = 90;', 'elements = 50;');
fid = fopen('TutorialShells50.m', 'w');
fwrite(fid, code50);
fclose(fid);

code70 = strrep(code, 'classdef TutorialShells < handle', 'classdef TutorialShells70 < handle');
code70 = strrep(code70, 'function obj = TutorialShells()', 'function obj = TutorialShells70()');
code70 = strrep(code70, 'elements = 90;', 'elements = 70;');
fid = fopen('TutorialShells70.m', 'w');
fwrite(fid, code70);
fclose(fid);


try
    disp('Running 50x50...');
    obj50 = TutorialShells50();
    
    disp('Running 70x70...');
    obj70 = TutorialShells70();
    
    w50 = obj50.wFun.fValues;
    coords50 = obj50.mesh.coord;
    
    % Find leading edge tip. The max deflection usually occurs there.
    % We just pick the node with the maximum absolute displacement.
    [maxW50, i50] = max(abs(w50));
    val50 = w50(i50);
    
    w70 = obj70.wFun.fValues;
    coords70 = obj70.mesh.coord;
    [maxW70, i70] = max(abs(w70));
    val70 = w70(i70);
    
    fprintf('50x50: Max W = %.15e at (%.4f, %.4f)\n', val50, coords50(i50,1), coords50(i50,2));
    fprintf('70x70: Max W = %.15e at (%.4f, %.4f)\n', val70, coords70(i70,1), coords70(i70,2));
    
    error = abs(val70 - val50) / abs(val70) * 100;
    fprintf('Relative Error = %.6f %%\n', error);
catch ME
    disp(ME.message);
end

exit;
