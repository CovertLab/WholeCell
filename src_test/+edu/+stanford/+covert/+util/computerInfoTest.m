classdef computerInfoTest < TestCase
    methods
        function this = computerInfoTest(name)
            this = this@TestCase(name);
        end
        
        function testComputerInformation(~)
            [userName,hostName,ipAddress] = ...
                edu.stanford.covert.util.computerInfo();
            assertFalse(isempty(userName));
            assertFalse(...
                isempty(regexp(hostName, '[a-zA-Z0-9_-]+', 'once')),...
                ['Host name: ' hostName]);
            assertFalse(...
                isempty(regexp(ipAddress, '\d{1,3}\.\d{1,3}\.\d{1,3}\.\d{1,3}', 'once')),...
                ['IP Address: ' ipAddress]);
        end
    end
end
