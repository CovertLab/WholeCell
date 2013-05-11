function tf = strlexne(a, b)
%STRLEXNE Lexicographic not equal to.
%
%   STRLEXNE(A, B) returns 1 if A is lexicographically not equal to B, and 0
%   otherwise.
%
%   This is a MATLAB version of the Perl `ne' operator.

%   See also NE.

%   Author:      Peter J. Acklam
%   Time-stamp:  2003-10-13 11:15:37 +0200
%   E-mail:      pjacklam@online.no
%   URL:         http://home.online.no/~pjacklam

   tf = strlexcmp(a, b) ~= 0;
