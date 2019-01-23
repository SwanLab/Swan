function [ phifunct ] = ini_phifunct(npnod,coordinates)

D = 2*0.5;
phifunct = -1/sqrt(D)*ones(npnod,1);

% phifunct = 0.6*ones(npnod,1);
% for i=1:npnod
%     x = coordinates(i,1); y = coordinates(i,2);
%     if (y<=0.2 || y>=0.8)
%         phifunct(i) = -0.4;
%     end
% end


end

