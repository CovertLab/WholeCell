% PDF Utility functions
%
% Author: Jonathan Karr
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 10/5/2011
classdef PdfUtil
    methods (Static = true)
        function rotatePages(fileName)
            if ispc
                path = 'lib\pdftk\bin\';
            else
                path = '';
            end
            
            [status, result] = system([path 'pdftk ' fileName ' cat 1N output tmp.pdf']);
            if status
                throw(MException('PdfUtil:pdftkError', result));
            elseif ~exist('tmp.pdf', 'file')
                throw(MException('PdfUtil:pdftkError', 'pdftk did not convert file'));
            end
            
            delete(fileName);
            movefile('tmp.pdf', fileName);
        end
    end
end