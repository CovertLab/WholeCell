function sout = ddeblank(s)
%DDEBLANK Double deblank. Strip both leading and trailing blanks.
%
%   DDEBLANK(S) removes leading and trailing blanks and null characters from
%   the string S.  A null character is one that has a value of 0.
%
%   See also DEBLANK, DEWHITE, DDEWHITE.

%   Author:      Peter J. Acklam
%   Time-stamp:  2003-10-13 11:13:07 +0200
%   E-mail:      pjacklam@online.no
%   URL:         http://home.online.no/~pjacklam

   error(nargchk(1, 1, nargin));
   if ~ischar(s)
      warning('Input must be a string (char array).');
   end

   if isempty(s)
      sout = s;
      return;
   end

   [r, c] = find( (s ~= ' ') & (s ~= 0) );
   if size(s, 1) == 1
      sout = s(min(c) : max(c));
   else
      sout = s(:, min(c) : max(c));
   end
