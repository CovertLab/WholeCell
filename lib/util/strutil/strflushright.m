function lines = strflushright(lines, offset)
%STRFLUSHRIGHT Flush right each line of text.
%
%   STRFLUSHRIGHT(LINES, OFFSET) where LINES is a cell array of lines and
%   OFFSET is a non-negative integer, does a right flush of the text in each
%   line so the last character in each line is at the specified offset.  If
%   OFFSET is omitted, 72 is used.

%   Author:      Peter J. Acklam
%   Time-stamp:  2003-10-13 11:13:53 +0200
%   E-mail:      pjacklam@online.no
%   URL:         http://home.online.no/~pjacklam

   error(nargchk(1, 2, nargin));
   if nargin < 2
      offset = 72;
   end

   spc = ' ';
   tab = sprintf('\t');
   fmt = sprintf('%%%ds', offset);
   lines = lines(:);

   for i = 1:length(lines)
      k = find(lines{i} ~= spc & lines{i} ~= tab);
      if isempty(k)
         lines{i} = '';
      else
         lines{i} = sprintf(fmt, lines{i}(k(1):k(end)));
      end
   end
