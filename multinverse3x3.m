function [ b,det ] = multinverse3x3( a )
% computation of the inverse of 3x3 matrix (simultaneously)

% d : traspose of the adjoint matrix
d(1,1,:) = a(2,2,:).*a(3,3,:)-a(2,3,:).*a(3,2,:);
d(1,2,:) = a(3,2,:).*a(1,3,:)-a(3,3,:).*a(1,2,:);
d(1,3,:) = a(1,2,:).*a(2,3,:)-a(1,3,:).*a(2,2,:);

d(2,1,:) = a(2,3,:).*a(3,1,:)-a(2,1,:).*a(3,3,:);
d(2,2,:) = a(3,3,:).*a(1,1,:)-a(3,1,:).*a(1,3,:);
d(2,3,:) = a(1,3,:).*a(2,1,:)-a(1,1,:).*a(2,3,:);

d(3,1,:) = a(2,1,:).*a(3,2,:)-a(2,2,:).*a(3,1,:);
d(3,2,:) = a(3,1,:).*a(1,2,:)-a(3,2,:).*a(1,1,:);
d(3,3,:) = a(1,1,:).*a(2,2,:)-a(1,2,:).*a(2,1,:);

% determinant : det(A) = tr(A*adj(A))/n
det = zeros(size(a,3),1);
for i=1:3
    for k=1:3
        det = det + squeeze(a(i,k,:).*d(k,i,:));
    end
end
det = det/3;

% inverse of the matrix 'a', b = 1/det(A) * adj(A)
for i=1:3
    for j=1:3
        b(i,j,:)=squeeze(d(i,j,:))./det;
    end
end
end


