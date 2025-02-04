clc,clear,close all
s.E          = 1;
s.nu         = 0.3;
s.meshType   = 'Square';
s.meshN      = 150;
s.holeType   = 'Ellipse';
s.nSteps     = [20,20];
s.pnorm      = 'Inf';
s.damageType = "Area";
PFH = PhaseFieldHomogenizer(s);
[mat,phi,holeParam] = PFH.computeHomogMaterial();
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

for i=1:9
    subplot(3,3,i)
    i1 = ceil(i/3);
    i2 = rem(i-1,3)+1;
    if ismember(i,[1,2,4,5,9])
    surf(holeParam{1},holeParam{2},squeeze(mat(i1,i2,:,:)))
    xlabel('m2')
    ylabel('m1')
    end
end
title('Ellipse homogenization')