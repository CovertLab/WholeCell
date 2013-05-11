% edu.stanford.covert.io.editMatFileHeader and
% edu.stanford.covert.io.readMatFileHeader test cases
%
% Author: Jonathan Karr
% Author: Jared Jacobs
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 9/14/2010
classdef editMatFileHeaderTest < TestCase
    methods
        function this = editMatFileHeaderTest(name)
            this = this@TestCase(name);
        end
        
        function testEditMatFileTimestamp(~)
            %generate mat file
            matFile = [tempname '.mat'];
            save(matFile, 'matFile', 'matFile');
            timestamp = edu.stanford.covert.io.readMatFileHeader(matFile);
            assertFalse(733899 == timestamp);

            %change timestamp
            edu.stanford.covert.io.editMatFileHeader(matFile, 733899, 'PCWIN', '5.0');
            [timestamp, platform, version] = ...
                edu.stanford.covert.io.readMatFileHeader(matFile);
            assertEqual(733899, timestamp);
            assertEqual('PCWIN', platform);
            assertEqual('5.0', version);

            %cleanup
            delete(matFile);
        end
    end
end
