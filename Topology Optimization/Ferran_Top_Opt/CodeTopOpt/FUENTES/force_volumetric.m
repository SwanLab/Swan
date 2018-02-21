function f = force_volumetric(coord,ndime,ftype)


switch ftype
    case {'ELASTIC'}

    if ndime == 2
        x = squeeze(coord(1,:,:));
        y = squeeze(coord(2,:,:));

        f(1,:,:) = 0*(x.^2).*y;
        f(2,:,:) = 0*(y.^2);  
    end


    if ndime == 3
        x = squeeze(coord(1,:,:));
        y = squeeze(coord(2,:,:));
        z = squeeze(coord(2,:,:));

        f(1,:,:) = 0;  
        f(2,:,:) = 0;
        f(3,:,:) = 0;
    end

    case {'THERMAL'}
        if ndime == 1
            x = squeeze(coord(1,:,:));
                      
            %f(1,:,:) = -(x.^2 );
            f(1,:,:) = ones(size(x));
        end

     if ndime == 2
        x = squeeze(coord(1,:,:));
        y = squeeze(coord(2,:,:));

        %f(1,:,:) = -(x.^2 - y.^2);
        %f(1,:,:) = x.*y;
        f(1,:,:) = 1*ones(size(x'));
       % f(1,:,:) = zeros(size(x));
        
    end


    if ndime == 3
        x = squeeze(coord(1,:,:));
        y = squeeze(coord(2,:,:));
        z = squeeze(coord(2,:,:));

        f(1,:,:) = 0;  
    end
 

end








end