s.E = 210;
s.nu = 0.3;
s.meshType  = 'Square';
s.meshN     = 100;
s.holeType  = "Ellipse";
s.nSteps     = [5 5];
s.damageType = "Area";
PFH = PhaseFieldHomogenizer(s);
[mat,phi] = PFH.computeHomogMaterial();



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








