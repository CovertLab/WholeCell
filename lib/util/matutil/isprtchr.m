function t = isprtchr(c)
%ISPRTCHR True for printable characters.
%
%   For a string C, ISPRTCHR(C) is 1 for printable characters and 0
%   otherwise.
%
%   See also: ISALNUM, ISALPHA, ISASCII, ISDIGIT, ISLOWER, ISUPPER,
%   ISXDIGIT.

%   Author:      Peter J. Acklam
%   Time-stamp:  2002-03-03 13:50:57 +0100
%   E-mail:      pjacklam@online.no
%   URL:         http://home.online.no/~pjacklam

   error(nargchk(1, 1, nargin));

   t = ischar(c) & (   ( (  32 <= c ) & ( c <= 126 ) ) ...
                     | ( ( 159 <= c ) & ( c <= 255 ) ) );

   % whether this is correct depends on the coding system used
%   if strncmp(computer, 'PC', 2)
%      t = t | ( ( 130 <= c ) & ( c <= 140 ) ) ...
%            | ( ( 145 <= c ) & ( c <= 156 ) );
%   end
