function [ sid ] = get_mac_adress

not_win = true;
switch computer('arch')
    case {'maci','maci64'}
        [~,mac_add] = system('ifconfig |grep ether | grep -o -E "([[:xdigit:]]{1,2}:){5}[[:xdigit:]]{1,2}"');
    case {'glnx86','glnxa64'}
        [~,mac_add] = system('ifconfig  | grep HWaddr | grep -o -E "([[:xdigit:]]{1,2}:){5}[[:xdigit:]]{1,2}"');
    case {'win32','win64'}
        not_win = false;
        sid = '';
        ni = java.net.NetworkInterface.getNetworkInterfaces;
        while ni.hasMoreElements
            addr = ni.nextElement.getHardwareAddress;
            if ~isempty(addr)
                sid = [sid, '.', sprintf('%.2X', typecast(addr, 'uint8'))];
            end
        end
    otherwise, error('Unknown architecture');
end

if(not_win)
    mac_add = regexprep(mac_add,'\r\n|\n|\r','.');
    sid = upper(strrep(mac_add(1:end-1),':',''));
end


end

