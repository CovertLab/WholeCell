function tf = strlexlt(a, b)
%STRLEXLT Lexicographic less than.
%
%   STRLEXLT(A, B) returns 1 if A is lexicographically less than B, and 0
%   otherwise.
%
%   This is a MATLAB version of the Perl `lt' operator.
%
%   See also LT.

%   Author:      Peter J. Acklam
%   Time-stamp:  2003-10-13 11:15:50 +0200
%   E-mail:      pjacklam@online.no
%   URL:         http://home.online.no/~pjacklam

   tf = strlexcmp(a, b) < 0;
