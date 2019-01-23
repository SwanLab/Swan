function F_obj_and_theta_vademecum
Path = '/media/Alex/Vademecum';
%isnap = 1;
for isnap = 1:2145
fid = fopen([Path,'/',num2str(isnap),'/executed_cases.txt'],'r');
C = textscan(fid,repmat(['%s%f'],1,7));
fclose(fid);  
columns = length(C{1});
fobj = C{4};
theta = C{6};
Volum = C{10};
%(1)
figure(1)
plot(fobj);
figure(2)
plot(theta)
figure(3)
plot(Volum)

end

%
%C = textscan(lastline, '%s%s%s%f%s%f%s%f%s%f%s%f%s%f%s%f%s%f%s%f%s%f%s%f%s%f%s%f');
%    k = 1;
%    for i = 2:2:28
%        Datos{k} = [C{i}];
%        k = k + 1;
%    end
%
%
%figure(1)
%  plot(A);
%figure(2)
%plot(theta);
%
%
%end
%
