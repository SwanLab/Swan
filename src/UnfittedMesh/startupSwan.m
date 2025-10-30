% Set up startup file
% 1. Change userpath to your Swan path
% 2. Copy this script and change the name to "startup.m"
% 3. Move the script to "~/Matlab/R202XX/toolbox/local"
% 4. In Matlab, go to "preferences" and in "Matlab/General" 
%    click on "Update Toolbox Path Cache"
%
% With this, everytime Matlab is opened, Swan folder and subfolders
% will be included. If Github branch is changed, just run commmand 
% "startup"
userpath = 'C:\Users\giova\Swan\Swan';
addpath(genpath(userpath));
rmpath(genpath('Old'))
cd(userpath);