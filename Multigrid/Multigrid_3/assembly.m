function [Ak,bk] = assembly(p)
%builds the assembly matrix and returns associated vector.  
  
%Here we require the entire vectors of p
%we only allow a kth triangle of our mesh 
Ak = sparse(zeros(3));
id = eye(3); c = [];
bk = zeros(3,1);
A = [1, p(1,1), p(1,2);...
     1, p(2,1), p(2,2);...
     1, p(3,1), p(3,2)];

for k = 1:3  
    b = id(:,k);
    d = A\b;
    c(:,end+1) = d;
end

AreaK = polyarea(p(1:3,1),p(1:3,2));

for i = 1:3
    for j = 1:3       
        Ak(i,j) = AreaK * (c(2,i)*c(2,j) + c(3,i)*c(3,j));   
    end 
    bk(i,1) = AreaK/3;
end

end


