function t = isalpha(c)
%ISALPHA True for letters.
%
%   For a string C, ISALPHA(C) is 1 for letters and 0 otherwise.
%
%   See also: ISALNUM, ISASCII, ISDIGIT, ISLOWER, ISPRTCHR, ISUPPER,
%   ISXDIGIT.

%   Author:      Peter J. Acklam
%   Time-stamp:  2002-03-03 13:51:09 +0100
%   E-mail:      pjacklam@online.no
%   URL:         http://home.online.no/~pjacklam

   error(nargchk(1, 1, nargin));

   t = isletter(c);
