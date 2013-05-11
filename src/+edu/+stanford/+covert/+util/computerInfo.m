function [userName, hostName, ipAddress] = computerInfo()
%COMPUTERINFO Retrieves information about the machine running the code.

    h = java.net.InetAddress.getLocalHost();
    hostName = char(h.getHostName());
    ipAddress = char(h.getHostAddress().toString());
    switch computer
        case {'GLNX86','GLNXA64','MACI','MACI64'}
            userName = getenv('USER');
        case {'PCWIN','PCWIN64'}
            userName = getenv('USERNAME');
        otherwise
            throw(MException('computerInfo:UnsupportedArchitecture', computer));
    end
end
