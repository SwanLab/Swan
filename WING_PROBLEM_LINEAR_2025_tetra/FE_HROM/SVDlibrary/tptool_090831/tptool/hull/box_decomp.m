function [W B] = box_decomp(U)
%BOX_DECOMP

% TODO: doc
% centering is recommended

% 2^n vertex!

M = [min(U); max(U)];
L = M(2,:) - M(1,:);

[n1 n2] = size(U);

%%% TODO: this is a hack
% if n2==2
% 	A = [linspace(0,1,n1)' linspace(1,0,n1)'];
% 	if norm(U*(U'*A)-A) < 1e-10
% 		% linear
% 		W = A;
% 		B = inv(U'*A);
% 		return
% 	end
% end

m = 2^n2;
B = zeros(m, n2);

% TODO: make it less ugly + efficient
r = 0:2:2*n2-2;
z = ones(1,n2);
z(1) = 0;
for i = 1:m
	z = nexti(z);
	B(i,:) = M(z+r);
end

W = zeros(n1,m);
for j = 1:n1
	q = [(M(2,:)-U(j,:))./L; (U(j,:)-M(1,:))./L];
	z = ones(1,n2);
	z(1) = 0;
	for i = 1:m
		z = nexti(z);
		W(j,i) = prod(q(z+r));
	end
end
%W = [U ones(size(U,1),1)]/[B ones(m,1)];

% next index
function z = nexti(z)
i = 1;
z(i) = z(i) + 1;
while z(i) > 2
	z(i) = 1;
	i = i + 1;
	z(i) = z(i) + 1;
end
