function t = ismatrix(x)
%ISMATRIX True for matrix input.
%
%   ISMATRIX(X) returns 1 if it's argument is a matrix and 0 otherwise.  An
%   array is considered a matrix if it has exactly two dimensions.

%   Author:      Peter J. Acklam
%   Time-stamp:  2002-03-03 13:51:01 +0100
%   E-mail:      pjacklam@online.no
%   URL:         http://home.online.no/~pjacklam

   t = ndims(x) == 2;
