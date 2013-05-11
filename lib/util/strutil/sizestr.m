function str = sizestr(array, before, between, after)
%SIZESTR Return size of array as a string.
%
%   STR = SIZESTR(ARRAY, BEFORE, BETWEEN, AFTER) returns the size of the array
%   ARRAY as a string, where BEFORE is the string written before the first size
%   value, BETWEEN is the string written between each size value, and AFTER is
%   the string written after the last size value.
%
%   Examples:
%
%      sizestr(rand(3,1,4), '[', 'x', ']')      returns '[3x1x4]'
%      sizestr(rand(3,1,4), '', '-by-', '')     returns '3-by-1-by-4'
%      sizestr(rand(3,1,4), '( ', ' by ', ' )') returns '( 3 by 1 by 4 )'

%   Author:      Peter J. Acklam
%   Time-stamp:  2003-10-13 11:13:45 +0200
%   E-mail:      pjacklam@online.no
%   URL:         http://home.online.no/~pjacklam

   error(nargchk(4, 4, nargin));

   siz = size(array);

   fmt = sprintf('%s%%d', between);
   str = sprintf(fmt, siz(2:end));
   str = sprintf('%s%d%s%s', before, siz(1), str, after);
