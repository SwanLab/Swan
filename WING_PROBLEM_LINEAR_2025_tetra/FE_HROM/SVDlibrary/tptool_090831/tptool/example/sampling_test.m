D = sampling_func(@(x) 1, [], []);
if D ~= 1
	error
end

D = sampling_func(@(x) x^2, [-1 1], 3);
if ~all(D == [1; 0; 1])
	error
end

D = sampling_func(@(x) x(1)+x(2), [-1 1; 0 1], [3 2]);
if ~all(all(D == [-1 0; 0 1; 1 2]))
	error
end

D = sampling_func(@(x) x(1)+x(2)+x(3), [-1 1; 0 1; 2 3], [2 2 2]);
Dok = zeros(2,2,2);
Dok(:,:,1) = [1 2; 3 4];
Dok(:,:,2) = [2 3; 4 5];
if ~all(all(all(D == Dok)))
	error
end

D = sampling_vec(@(x) [x(1); x(2); 3], [-1 1; 0 1], [2 2]);
Dok = zeros(2,2,3);
Dok(1, 1, :) = [-1; 0; 3];
Dok(1, 2, :) = [-1; 1; 3];
Dok(2, 1, :) = [ 1; 0; 3];
Dok(2, 2, :) = [ 1; 1; 3];
if ~all(all(D == Dok))
	error
end

lpv = {@(x)1 @(x)x(2); @(x)0 @(x)x(1)};
dep = zeros(2,2,2);
dep(1,2,:)=[0 1];
dep(2,2,:)=[1 0];
data = sampling_lpv(lpv, dep, [-1 1; 4 5], [5 7]);
if data{1,1} ~= 1 || data{2,1} ~= 0
	error
end
