function v = regularize3(x,V_vect,dx,dy,dz)
% REGULARIZATION OF THE ADVECTION VELOCITY
% Here we regularize the velocity field by considering the H1 product instead of the L2 one.

%   k1 = e2/(dx^2) ; % Coefficients of the regularization matrix.
%   k2 = e2/(dy^2) ;
%   k3 = e2/(dz^2) ;

% !! PATCH !!
load(fullfile(pwd,'Allaire_ShapeOpt','conversion'));

for n = 1:length(x)
    phi(b1(n,1),b1(n,2),b1(n,3)) = x(n);
end

for n = 1:length(V_vect)
    V(b1(n,1),b1(n,2),b1(n,3)) = V_vect(n);
end
% figure, surf(-permute(V(:,ceil(size(V,2)/2),:),[1 3 2])), title('Raw V')
% figure('NumberTitle', 'off', 'Name', 'FEM-MAT-OO - Raw V')
% subplot(2,2,1), surf(-V(:,:,2)), title('V - Root')
% subplot(2,2,2), surf(-V(:,:,end-1)), title('V - Tip')
% subplot(2,2,3), surf(permute(-V(ceil(size(V,1)/2),:,:),[2 3 1])), title('V - XY')
% subplot(2,2,4), surf(permute(-V(:,ceil(size(V,2)/2),:),[1 3 2])), title('V - XZ')

% Now we calculate the surface Dirac function.
epsperim =min( min(dx,dy),dz)/20 ;
sx = phi./sqrt(phi.^2+epsperim^2) ;

sxn = shift3n('n',sx,'ng') ;
sxs = shift3n('s',sx,'ng') ;
sxe = shift3n('e',sx,'ng') ;
sxw = shift3n('w',sx,'ng') ;
sxu = shift3n('u',sx,'ng') ;
sxd = shift3n('d',sx,'ng') ;

%We now calculate d(phi)/dx and d(phi)/dy:
dsxx = (sxw-sxe)/(2*dx) ;
dsxy = (sxn-sxs)/(2*dy) ;
dsxz = (sxu-sxd)/(2*dz) ;

delta = 0.5*sqrt(dsxx.^2+dsxy.^2+dsxz.^2); % Surface Dirac function

b = -V.*delta; % The right hand side is defined on the boundary so we use the delta function.

b_vect(A1(:,:,:)) = b(:,:,:);
b_vect = b_vect';


filter = Filter_PDE_Density('Cantileverbeam_Hexahedra_Bilinear_Structured','MACRO');
filter.preProcess;
filter.updateEpsilon(0.03);

v = filter.getP1fromP1(b_vect);
for n = 1:length(v)
    V(b1(n,1),b1(n,2)) = v(n);
end

% figure, surf(-permute(V(:,ceil(size(V,2)/2),:),[1 3 2])), title('Regularized V')

% figure('NumberTitle', 'off', 'Name', 'FEM-MAT-OO - Regularized V')
% subplot(2,2,1), surf(-V(:,:,2)), title('V - Root')
% subplot(2,2,2), surf(-V(:,:,end-1)), title('V - Tip')
% subplot(2,2,3), surf(permute(-V(ceil(size(V,1)/2),:,:),[2 3 1])), title('V - XY')
% subplot(2,2,4), surf(permute(-V(:,ceil(size(V,2)/2),:),[1 3 2])), title('V - XZ')

end