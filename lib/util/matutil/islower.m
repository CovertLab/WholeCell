function t = islower(c)
%ISLOWER True for lowercase letters.
%
%   For a string C, ISLOWER(C) is 1 for lowercase letters and 0 otherwise.
%
%   See also: ISALNUM, ISALPHA, ISASCII, ISDIGIT, ISPRTCHR, ISUPPER,
%   ISXDIGIT.

%   Author:      Peter J. Acklam
%   Time-stamp:  2002-03-03 13:51:03 +0100
%   E-mail:      pjacklam@online.no
%   URL:         http://home.online.no/~pjacklam

   error(nargchk(1, 1, nargin));

   t = isletter(c) & (c == lower(c));
