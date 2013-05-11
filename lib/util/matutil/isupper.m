function t = isupper(c)
%ISUPPER  True for uppercase letters.
%
%   For a string C, ISUPPER(C) is 1 for uppercase letters and 0
%   otherwise.
%
%   See also: ISALNUM, ISALPHA, ISASCII, ISDIGIT, ISLOWER, ISPRTCHR,
%   ISXDIGIT.

%   Author:      Peter J. Acklam
%   Time-stamp:  2002-03-03 13:50:48 +0100
%   E-mail:      pjacklam@online.no
%   URL:         http://home.online.no/~pjacklam

   error(nargchk(1, 1, nargin));

   t = isletter(c) & (c == upper(c));
