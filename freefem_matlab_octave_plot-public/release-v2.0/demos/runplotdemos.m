clear all;
fprintf('capacitor_2d.edp\n');
system('FreeFem++ -ns -nw -v 0 capacitor_2d.edp');
fprintf('interpolate_complex.edp\n');
system('FreeFem++ -ns -nw -v 0 interpolate_complex.edp');
fprintf('capacitor_3d.edp\n');
system('FreeFem++ -ns -nw -v 0 capacitor_3d.edp');
fprintf('region.edp\n');
system('FreeFem++ -ns -nw -v 0 regions.edp');
fprintf('periodic_bc.edp\n');
system('FreeFem++ -ns -nw -v 0 periodic_bc.edp');
fprintf('heat_transfer_p1b.edp\n');
system('FreeFem++ -ns -nw -v 0 heat_transfer_p1b.edp');
fprintf('convective_rolls.edp\n');
system('FreeFem++ -ns -nw -v 0 convective_rolls.edp');
fprintf('demo_Lshape.edp\n');
system('FreeFem++ -ns -nw -v 0 demo_Lshape.edp');
fprintf('vectored_fespace.edp\n');
system('FreeFem++ -ns -nw -v 0 vectored_fespace.edp');
capacitor_2d
pause(5);
close all;
gui
close all;
bdcoloring
pause(5);
close all;
interpolate
pause(5);
close all;
capacitor_3d
pause(5);
close all;
region
pause(5);
close all;
periodic_bc
pause(5);
close all;
heat_transfer
pause(5);
close all;
convective_rolls
pause(5);
close all;
demo_Lshape
pause(5);
close all;
vectored_fespace
pause(5);
close all;

