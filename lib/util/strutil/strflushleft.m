function lines = strflushleft(lines, offset)
%STRFLUSHLEFT Flush left each line of text.
%
%   STRFLUSHLEFT(LINES, OFFSET) where LINES is a cell array of lines and OFFSET
%   is a non-negative integer, does a left flush of the text in each line and
%   indents it to the specified offset.

%   Author:      Peter J. Acklam
%   Time-stamp:  2003-10-13 11:12:55 +0200
%   E-mail:      pjacklam@online.no
%   URL:         http://home.online.no/~pjacklam

   error(nargchk(1, 2, nargin));
   if nargin < 2
      offset = 0;
   end

   spc = ' ';
   tab = sprintf('\t');
   pad = spc(:,ones(offset,1));
   lines = lines(:);

   for i = 1:length(lines)
      k = find(lines{i} ~= spc & lines{i} ~= tab);
      if isempty(k)
         lines{i} = '';
      else
         lines{i} = [pad lines{i}(k(1):k(end))];
      end
   end
