% Set up startup file
% 1. Change userpath to your Swan path
% 2. Copy this script in "~/Matlab/R202XX/toolbox/local"
% 3. Change file name to "startup.m"
% 4. In Matlab, go to "preferences" and in "Matlab/General" 
%    click on "Update Toolbox Path Cache"
%
% With this, everytime Matlab is opened, Swan folder and subfolders
% will be included. If Github branch is changed, just run commmand 
% "addpath(genpath(userpath))"
userpath = 'C:\Users\XXXX\Documents\GitHub\Swan';
addpath(genpath(userpath));
cd(userpath);