%Path = '/home/aferrer/Dropbox/LisboaCongres/';
%newmicros = [23,4758,4523,31,6115]';
Path = '/home/aferrer/Desktop/Micros/';
newmicros = randi([1 8000],100,1);
iter = 1;
post_micro = 'images';
postprocess_micro(Path,newmicros,iter)
postprocess_micros(Path,newmicros,post_micro)