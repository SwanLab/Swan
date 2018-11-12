function symetries_prueba
n_cuadrantes = 4;
phi_desv = pi/12;
theta_desv = pi/12;

phi_cuadr = [-pi/4 -3*pi/4 3*pi/4 pi/4]';
phi_minus = -phi_desv*ones(4,1);
phi_plus = phi_desv*ones(4,1);

phi_all = [phi_cuadr+phi_plus,phi_cuadr+phi_minus]';
phi_all = phi_all(:);
theta_all = pi/2-theta_desv*ones(size(phi_all));
strain = phitheta2strain(phi_all,theta_all);
A = zeros(3,3,size(strain,1));
for i = 1:size(strain,1)
A(:,:,i) = strain(i,:)'*strain(i,:);
end

A(:,:,[1])
A(:,:,[8])

Ch1 = read_Ch('/home/aferrer/Documents/Doctorat/Tesi/MicroEscalaDT/Cases4FindingSymmetries/1/Ch.txt')
Ch7 = read_Ch('/home/aferrer/Documents/Doctorat/Tesi/MicroEscalaDT/Cases4FindingSymmetries/8/Ch.txt')




end



function Ch_mat = read_Ch(Path)
fid = fopen(Path,'r');
C = textscan(fid,repmat('%s%f',1,11));
Ch = C{4:2:14};
C11 = C{4};
C12 = C{6};
C13 = C{8};
C22 = C{10};
C23 = C{12};
C33 = C{14};
Ch_mat = [C11(end) C12(end) C13(end); C12(end) C22(end) C23(end); C13(end) C23(end) C33(end)];
fclose(fid);

end

function strain = phitheta2strain(phi,theta)

stra_theta = theta;
stra_phi = phi;
stra_x = sin(stra_theta).*cos(stra_phi);
stra_y = sin(stra_theta).*sin(stra_phi);
stra_xy = cos(stra_theta);
strain = [stra_x stra_y stra_xy];

end
