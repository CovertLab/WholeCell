function t = isalnum(c)
%ISALNUM True for alphanumerics (letters and digits).
%
%   For a string C, ISALNUM(C) is 1 for alphanumerics (letters and
%   digits) and 0 otherwise.
%
%   See also: ISALPHA, ISASCII, ISDIGIT, ISLOWER, ISPRTCHR, ISUPPER,
%   ISXDIGIT.

%   Author:      Peter J. Acklam
%   Time-stamp:  2002-03-03 13:51:11 +0100
%   E-mail:      pjacklam@online.no
%   URL:         http://home.online.no/~pjacklam

   error(nargchk(1, 1, nargin));

   t = isletter(c) | isdigit(c);
