function tf = strlexge(a, b)
%STRLEXGE Lexicographic greater than or equal to.
%
%   STRLEXGE(A, B) returns 1 if A is lexicographically greater than or equal to
%   B, and 0 otherwise.
%
%   This is a MATLAB version of the Perl `ge' operator.
%
%   See also GE.

%   Author:      Peter J. Acklam
%   Time-stamp:  2003-10-13 11:16:04 +0200
%   E-mail:      pjacklam@online.no
%   URL:         http://home.online.no/~pjacklam

   tf = strlexcmp(a, b) >= 0;
