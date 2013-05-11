function tf = strlexle(a, b)
%STRLEXLE Lexicographic less than or equal to.
%
%   STRLEXLE(A, B) returns 1 if A is lexicographically less than or equal to B,
%   and 0 otherwise.
%
%   This is a MATLAB version of the Perl `le' operator.
%
%   See also LE.

%   Author:      Peter J. Acklam
%   Time-stamp:  2003-10-13 11:15:58 +0200
%   E-mail:      pjacklam@online.no
%   URL:         http://home.online.no/~pjacklam

   tf = strlexcmp(a, b) <= 0;
