function [cartd,djacb,injacb] = cal_cartd(igaus,posgp,element,ndime,nnode,nelem,coordinatesn,coordinatesa,...
    ptype)

etype = element.type;

neres = 0; %hessian
dirichlet_data = zeros(nnode,nelem);
coordn = zeros(nnode,ndime,nelem);
%coorda = zeros(nnode,ndime,nelem);
jacob = zeros(ndime,ndime,nelem);
injacb= zeros(ndime,ndime,nelem);
cartd = zeros(ndime,nnode,nelem);



[shape,deriv,heslo] = shape_deriv_functions(igaus,posgp,ptype,etype,nnode,neres);
% for i=1:nnode
%     dirichlet_data(i,:)= element.conectivities(:,i);
%     for idime=1:ndime
%         coordn(i,idime,:)= coordinatesn(dirichlet_data(i,:),idime);
%         %coorda(i,idime,:)= coordinatesa(dirichlet_data(i,:),idime);
%     end
% end

for idime = 1:ndime
    a = coordinatesn(:,idime);
    coordn(:,idime,:) = a(permute(element.conectivities',[1,3,2]));
end



elcod = coordn;
% for idime = 1:ndime
%    for jdime = 1:ndime
%        for inode = 1:nnode
%            jacob(idime,jdime,:) = jacob(idime,jdime,:) + deriv(idime,inode,:).*elcod(inode,jdime,:);
%        end
%    end
% end

for i_dime=1:ndime
    jacob(i_dime,:,:) = sum(repmat(permute(deriv(i_dime,:,:),[2,1,3]),1,ndime,nelem) .* elcod,1);
end


% inverse of jacobian matrix (inverse:injacb, det=djacb
switch ptype
    case '1D'
        injacb = 1./jacob;
        djacb = squeeze(jacob);
    case '2D'
        [injacb,djacb] = multinverse2x2(jacob);
    case '3D'
        [injacb,djacb] = multinverse3x3(jacob);
end

% cartesian derivates
% for idime = 1:ndime
%     for inode = 1:nnode
%         for jdime = 1:ndime
%             cartd(idime,inode,:) = cartd(idime,inode,:)+injacb(idime,jdime,:).*deriv(jdime,inode,:);
%         end
%     end
% end


for i_dime=1:ndime
    cartd(i_dime,:,:) = sum(repmat(permute(injacb(i_dime,:,:),[2,1,3]),1,nnode,1) .* repmat(deriv,1,1,nelem),1);
end



end

