function t = isascii(c)
%ISASCII True for decimal digits.
%
%   For a string C, ISASCII(C) is 1 for any ASCII character code between 0
%   and 127, inclusive, and 0 otherwise.
%
%   See also: ISALNUM, ISALPHA, ISDIGIT, ISLOWER, ISPRTCHR, ISUPPER,
%   ISXDIGIT.

%   Author:      Peter J. Acklam
%   Time-stamp:  2002-03-03 13:51:07 +0100
%   E-mail:      pjacklam@online.no
%   URL:         http://home.online.no/~pjacklam

   error(nargchk(1, 1, nargin));

   t = ischar(c) & ( 0 <= c ) & ( c <= 127 );
