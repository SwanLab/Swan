% RunNullSimulations

% s.gJFlow = 0.7;
% s.Vf     = 0.2;
% Sim1     = ThreeDimCantileverDensity(s);
% 
% s.gJFlow = 2;
% s.Vf     = 0.4;
% Sim2     = ThreeDimCantileverDensity(s);








% MBB Beam 40% Density: gJ=1 and gJ=2
% MBBBeamDensity(1);
% MBBBeamDensity(2);

% 2D Cantilever Beam LevelSet: gJ=0.2,1,2.5
% TwoDimCantilever(0.2);
% TwoDimCantilever(1);
% TwoDimCantilever(2.5);

% 3D Cantilever Beam LevelSet: gJ=0.5/1/2 (V0d4) + gJ=0.7 (V0d2)
% CREATE FOLDER
% ThreeDimCantilever(0.5,0.4);
% ThreeDimCantilever(1,0.4);
% ThreeDimCantilever(2,0.4);
ThreeDimCantilever(0.7,0.2);

% Gripper LevelSet: gJ=0.5, gJ=0.05
% CREATE FOLDER
% Gripper(0.05);
% Gripper(0.5);









% MBB Beam 40% LevelSet: gJ= 0.5/1/2;
% MBBBeam(0.5,0.015);
% MBBBeam(1,0.03);
% MBBBeam(2,0.1);


% MultiLoad Bridge LevelSet: gJ=10/1Load + gJ=10/3Loads + gJ=10/9Loads;
% PENDING
