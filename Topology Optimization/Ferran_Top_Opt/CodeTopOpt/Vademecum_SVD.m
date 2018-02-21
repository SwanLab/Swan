function [phi,theta,level_set] = Vademecum_SVD()
pint = 0; nlevel = 7; % 0 == no pintar, 1 == pintar
[~,phi,theta,numsnapshots,phi_matrix,theta_matrix,number_matrix,theta_step] = vademecum_spehere_mesh(nlevel,pint);
coordinates = read_coordinates;
read = 0;

if read
level_set = read_level_set;
else
a = load ('/home/aferrer/Documents/Doctorat/Tesi/MicroStructures_Vademecum/level_set.mat','level_set');
level_set = a.level_set;
end

%One level set example
tri = delaunay(coordinates(:,1),coordinates(:,2)); 
h = trisurf(tri,coordinates(:,1),coordinates(:,2),level_set(:,5));


indexb = ((theta >= -1e-3) & (theta <= 1e-3)) | ((theta >= pi/2-1e-3) & (theta <= pi/2+ 1e-3));
[~,orded_index] = sort(phi(indexb),'descend');
vec_index = 1:size(level_set,2);
index = vec_index(indexb);
orded2 = index(orded_index);

% for i = 1:length(orded2)
% tri = delaunay(coordinates(:,1),coordinates(:,2)); 
% h = trisurf(tri,coordinates(:,1),coordinates(:,2),level_set(:,orded2(i)));   
% axis([0 1 0 1 -4 4])
% print(['/home/aferrer/Documents/Doctorat/Tesi/Rastro/figures/GifLevelSetVademecum/iteracio_',num2str(i)],'-deps')
% end


%All snapshots
figure(1)
err_trunc = compute_trunc_error(level_set);
plot((err_trunc),'b-+');
xlabel('Bases')
ylabel('Truncation error')
print(['/home/aferrer/Documents/Doctorat/Tesi/Rastro/figures/TruncErrorAllVademecum'],'-depsc')

figure(2)
err_trunc = compute_trunc_error(level_set(:,orded2));
plot((err_trunc),'b-+');
xlabel('Bases')
ylabel('Truncation error')
print(['/home/aferrer/Documents/Doctorat/Tesi/Rastro/figures/TruncErrorNoThetaVademecum'],'-depsc')

end

function coordinates = read_coordinates
path = '/home/aferrer/Documents/Doctorat/Tesi/MicroStructures_Vademecum/Vademecum_20_penalty_0_999_1gpFULL/';
[~,file] = unix(['ls ',path,num2str(1),'/*.msh']);    
fid = fopen(file(1:end-1),'r');
inode = 1;
for iline = 1:3287
    tline = fgets(fid);
    if iline > 6
        carac = str2num(tline);
        coordinates(inode,:) = carac(2:3);
        inode = inode +1;
    end
        
end
end


function err_trunc = compute_trunc_error(data)
[S,D,V] = svd(data);
sigm = diag(D);
c1 = sum(sigm.^2);
R = length(sigm);
err_trunc = zeros(R,1);
%err_trunc2 = zeros(R,1);
for i = 1: R
err_trunc(i) = sqrt((c1 - sum(sigm(1:i).^2))/c1);
%err_trunc2(i) = norm(level_set - S(:,1:i)*D(1:i,1:i)*V(:,1:i)','fro')/sqrt(c1);
end
end