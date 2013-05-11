function t = isxdigit(c)
%ISXDIGIT True for hexadecimal digits.
%
%   For a string C, ISXDIGIT(C) is 1 for hexadecimal digits and 0
%   otherwise.
%
%   (In ASCII: characters [0-9], [A-F] or [a-f]).
%
%   See also: ISALNUM, ISALPHA, ISASCII, ISDIGIT, ISLOWER, ISPRTCHR,
%   ISUPPER.

%   Author:      Peter J. Acklam
%   Time-stamp:  2002-03-03 13:50:45 +0100
%   E-mail:      pjacklam@online.no
%   URL:         http://home.online.no/~pjacklam

   error(nargchk(1, 1, nargin));

   t = ischar(c) & (   ( ( '0' <= c ) & ( c <= '9' ) ) ...
                     | ( ( 'A' <= c ) & ( c <= 'F' ) ) ...
                     | ( ( 'a' <= c ) & ( c <= 'f' ) ) );
