function [ b,det ] = multinverse2x2( a )
% computation of the inverse of 2x2 matrix (simultaneously)

% d : traspose of the adjoint matrix
d(1,1,:) =  a(2,2,:);
d(1,2,:) = -a(1,2,:);
d(2,1,:) = -a(2,1,:);
d(2,2,:) =  a(1,1,:);
% determinant
det = squeeze(a(1,1,:).*a(2,2,:)-a(1,2,:).*a(2,1,:));

for i=1:2
    for j=1:2
        b(i,j,:) = squeeze(d(i,j,:))./det;
    end
end
end

