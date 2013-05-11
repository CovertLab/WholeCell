function tf = strlexeq(a, b)
%STRLEXEQ Lexicographic equal to.
%
%   STRLEXEQ(A, B) returns 1 if A is lexicographically equal to B, and 0
%   otherwise.
%
%   This is a MATLAB version of the Perl `eq' operator.
%
%   See also EQ.

%   Author:      Peter J. Acklam
%   Time-stamp:  2003-10-13 11:16:01 +0200
%   E-mail:      pjacklam@online.no
%   URL:         http://home.online.no/~pjacklam

   tf = strlexcmp(a, b) == 0;
