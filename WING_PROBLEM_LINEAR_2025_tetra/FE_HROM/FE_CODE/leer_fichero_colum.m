

% It reads a graph file (two columns)
%function [DATOS]=leer_fichero_colum(name,ncolum,nbasura,nrow,nbasura_After,abrir,fid,MoreGP,ngaus)
%                                               %%%%%%%%%%%%%%%%%%%%%
%                                               Not required
function [DATOS]=leer_fichero_colum(name,ncolum,nbasura,nrow,nbasura_After,abrir,fid,MoreGP,ngaus)

% Defaults
if nargin == 0
    name = '/home/joaquin/Desktop/aaa'
    ncolum = 2
    nbasura = 2;   
    nrow    = -1;   
    nbasura_After = 0;
    abrir = 1;
    fid = 3;
    MoreGP = 0; % More than one gauss point (for reading results, .res file)
    ngaus = 1;
elseif nargin == 2
    nbasura = 0;   
    nrow    = -1;  
    nbasura_After = 0;
    abrir = 1;
    fid = -1;
    MoreGP = 0; % More than one gauss point (for reading results, .res file)
    ngaus = 1;
elseif nargin == 3
    nrow    = -1;    
    nbasura_After = 0;
    abrir = 1;
    fid = -1;
    MoreGP = 0; 
    ngaus = 1;
elseif nargin ==4  
    
    nbasura_After= 0;
    abrir = 1;
    fid = -1;
    MoreGP = 0; 
    ngaus = 1;
elseif nargin ==5 
    abrir = 1;
    fid = -1;
    MoreGP = 0; 
    ngaus = 1;    
elseif nargin == 7
    MoreGP = 0; 
    ngaus = 1;    
end 

if abrir == 1
    fid=fopen([name],'r');
else
    % It's an input
end
if fid == -1
    disp('There s no file')
    DATOS = 0;
else
    for i = 1:nbasura 
        linea_before=fgets(fid)
    end
    
    if nrow == -1 % We don't know how many lines have to be read 
        AAAA=fscanf(fid,'%f');        
        nrow = length(AAAA)/ncolum;
        DATOS = length(AAAA);
        for irow = 1:nrow
            for icol = 1:ncolum
                DATOS(irow,icol) = AAAA((irow-1)*ncolum+icol);
            end
        end
    else
        if MoreGP == 0
        DATOS = fscanf(fid,'%f',[  ncolum nrow] );
        %AAAA=fscanf(fid,'%f',ncolum*nrow);
        DATOS = DATOS';
        else 
            % Another data structure . Something like this: 
            %  1  a1_1    b1_1   c1_1 ...
            %     a1_2    b1_2   c1_2 ...
            %     . . .. .......
            %     a1_ng   b1_ng  c1_ng ...
            %
            %  2  a2_1    b2_1   c2_1 
            %      ....................
            %
            %
            %
            % %%%%%%
            
            % Initialize 
            DATOS = zeros(nrow*ngaus,ncolum);
            
            for irow = 1:nrow                              
                data1 = fscanf(fid,'%f \n',[1 ncolum]);
                data2 = fscanf(fid,'%f',[ncolum-1,ngaus-1]);             
                indI = (irow-1)*ngaus+1;
                indF = ngaus*irow;
                DATOS(indI,:) = data1;
                DATOS(indI+1:indF,2:end)=data2';
                DATOS(indI+1:indF,1)    = irow;
            end
            
        end
    end
    
    
    
    % Txt=fscanf(fid,'%s',WordAfter)
    for i = 1:nbasura_After+1 
        linea_After=fgets(fid);
    end
    
    

    
end

if abrir== 1
      sta = fclose(fid);
 
end

