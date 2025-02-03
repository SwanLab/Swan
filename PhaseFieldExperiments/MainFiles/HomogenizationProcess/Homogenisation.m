clc,clear,close all
s.E          = 210;
s.nu         = 0.3;
s.meshType   = 'Hexagon';
s.meshN      = 50;
s.holeType   = 'SmoothHexagon';
s.nSteps     = [20];
s.pnorm      = 'Inf';
s.damageType = "Area";
PFH = PhaseFieldHomogenizer(s);
[mat,phi] = PFH.computeHomogMaterial();
%save("RectangleMicroDamagePerimeter_min0,1","mat","phi")

%%%% EL RESULTAT S'HA DE DIVIDIR PEL VOLUM %%%%%%%%%

% [mat,phi] = PFH.computeHomogMaterial("Circle","Perimeter",100);
% save('CircleMicroDamagePerimeter','mat','phi')
% [mat,phi] = PFH.computeHomogMaterial("Circle","Area",100);
% save('CircleMicroDamageArea','mat','phi')
% 
% [mat,phi]  = PFH.computeHomogMaterial("Square","Perimeter",100);
% save('SquareMicroDamagePerimeter','mat','phi')
% [mat,phi]  = PFH.computeHomogMaterial("Square","Area",100);
% save('SquareMicroDamageArea','mat','phi')

% [mat,phi] = PFH.computeIsotropicMaterial("AT1",100);
% save('IsoMicroDamage','mat','phi')
% [~,AT2] = PFH.computeIsotropicMaterial("AT2",30);








