function t = isrowvector(x)
%ISROWVECTOR True for row vector input.
%
%   ISROWVECTOR(X) returns 1 if X is a row vector and 0 otherwise.
%
%   An array is considered a row vector if the length along each dimension,
%   except possibly the second, is never larger than one.  That is,
%   SIZE(X, DIM) <= 1 for all DIM except possibly when DIM = 2.

%   Author:      Peter J. Acklam
%   Time-stamp:  2002-03-03 13:50:55 +0100
%   E-mail:      pjacklam@online.no
%   URL:         http://home.online.no/~pjacklam

   t = ndims(x) == 2 & size(x, 1) == 1;
