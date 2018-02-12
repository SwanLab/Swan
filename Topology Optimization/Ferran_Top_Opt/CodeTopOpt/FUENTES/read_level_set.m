function level_set = read_level_set
%path = '/home/aferrer/Documents/Doctorat/Tesi/MicroStructures_Vademecum/Vademecum1/';
path = '/home/aferrer/Documents/Doctorat/Tesi/MicroStructures_Vademecum/Vademecum_20_penalty_0_999_1gpFULL/';
%path = /media/aferrer/Alex/Vademecum_20_penalty_0_999_1gp/';

[~,nsnap] = unix(['ls ',path,' | wc -l']);
nsnap = str2num(nsnap) -1;
for isnap = 1:nsnap
[~,file] = unix(['ls ',path,num2str(isnap),'/*.res']);    
fid = fopen(file(1:end-1),'r');

inode = 1;
for iline = 1:3294
    tline = fgets(fid);
    if iline > 13
        
        carac = str2num(tline);
       
        if isempty(carac)
            tline
        end
        level_set(inode,isnap) = carac(2);
        inode = inode +1;
    end
  
end
isnap
end

end