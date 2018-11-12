function [ working_path,results_path ] = switch_runner( Runner )

% Runner = char(java.net.InetAddress.getLocalHost.getHostName);
switch Runner
    case 'Ferran' %'FERRAN-PC'
        working_path = 'E:\Dropbox\Ferran\CodeTopOpt';
        results_path = 'E:\TFM\23 - Implementation new algorithms';
        
    case 'FerranLaptop' %'FERRAN-PORT'
        working_path = 'D:\Dropbox\Ferran\CodeTopOpt';
        results_path = 'D:\TFM\TestingMMA';
        
    case 'AlexLaptop'
        working_path = '/home/alex/Dropbox/Ferran/CodeTopOpt';
        results_path = '/home/alex/Desktop/Results/Results';
        
    case 'AlexOffice'
        working_path = '/home/aferrer/Dropbox/Ferran/CodeTopOpt';
        results_path = '/home/aferrer/Desktop/ResultsMicro/Results';
        
    otherwise
        error('Runner not detected.');
end


end

