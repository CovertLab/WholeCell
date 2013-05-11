function sout = dewhite(s)
%DEWHITE Dewhite. Strip trailing whitespace.
%
%   DEWHITE(S) removes leading and trailing white space and any null characters
%   from the string S.  A null character is one that has an absolute value of
%   0.
%
%   See also DDEWHITE, DEBLANK, DDEBLANK.

%   Author:      Peter J. Acklam
%   Time-stamp:  2003-10-13 11:12:52 +0200
%   E-mail:      pjacklam@online.no
%   URL:         http://home.online.no/~pjacklam

   error(nargchk(1, 1, nargin));
   if ~ischar(s)
      error( 'Input must be a string (char array).' );
   end

   if isempty(s)
      sout = s;
      return;
   end

   [r, c] = find(~isspace(s));
   if size(s, 1) == 1
      sout = s(1:max(c));
   else
      sout = s(:,1:max(c));
   end
