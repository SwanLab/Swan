clc
clear
close all

load('MS_meu.mat')
a_meu = a;
coord_meu = coord;
fV_meu = fV;
masters_meu = masters;
mR_meu = mR;
MS_meu = MS;
obj_meu = obj;
p1f_meu = p1f;
s_meu = s;
slaves_meu = slaves;


load('MS_seu.mat')
a_seu = a;
coord_seu = coord;
fV_seu = fV;
masters_seu = masters;
mR_seu = mR;
MS_seu = MS;
obj_seu = obj;
p1f_seu = p1f;
s_seu = s;
slaves_seu = slaves;

clear('a', 'coord', 'fV', 'masters', 'mR', 'MS', 'obj', 'p1f', 's', 'slaves')


if(isequal(a_seu,a_meu))
    a_equal = 'True';
else
    a_equal = 'False';
end

if(isequal(fV_seu,fV_meu))
    fV_equal = 'True';
else
    fV_equal = 'False';
end













