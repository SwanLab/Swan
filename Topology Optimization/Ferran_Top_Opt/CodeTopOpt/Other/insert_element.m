function [ result ] = insert_element( data,pos,x )

c = false(1,length(data)+length(pos));
c(pos) = true;
% result = nan(size(c));
result(~c) = data;
result(c) = x;

end

