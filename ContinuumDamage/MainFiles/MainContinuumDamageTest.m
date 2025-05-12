clc;clear;close all

%load('TestForceTraction1Elem.mat')
%load('TestDisplacementTraction.mat')
%cParams.mesh.name = 'CD_Mesh';

% type = '1Dtrac';
% 
% switch type
%     case 'SEMtrac'
%         cParams.mesh.name = 'PF_SENtraction0_0025';
%         cParams.bc.bcType = 'SEMtraction'; 
%     case 'SEMmix'
%         cParams.mesh.name = 'PF_SENmixed0_0025';
%         cParams.bc.bcType = 'SEMmixed'; 
%     case 'SEMshear'
%         cParams.mesh.name = 'PF_SENshear0_0025';
%         cParams.bc.bcType = 'SEMshear'; 
%     case '1Dtrac'
%         cParams.mesh.meshLength = 1;
%         cParams.mesh.meshWidth = 1;
%         cParams.mesh.meshN = 10;
%         cParams.mesh.meshM = 10;
%         cParams.bc.bcType = 'displacementTraction';
% end









% plotClass = ContinuumDamagePlotter(data);
% 
% plotClass.plotDisplacementField(); %FALTA POSAR LA H COM A INPUT
% plotClass.plotDamagesField();
% 
% disp = 'disp';
% dmg = 'max damage';
% frce = 'force';
% engy = 'total energy';
% rVar = 'max r';
% qVar = 'max q';
% mat = 'material';
% 
% plotClass.plotSelector (disp,dmg,'Damage - Displacement');
% plotClass.plotSelector (disp,frce,'Force - Displacement');
% plotClass.plotSelector (disp,engy,'Energy - Displacement');
% plotClass.plotSelector (rVar,qVar,'q - r');
% plotClass.plotSelector (rVar,dmg,'Damage-r');
% plotClass.plotSelector (rVar,frce,'Force-r');
