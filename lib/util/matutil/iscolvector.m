function t = iscolvector(x)
%ISCOLVECTOR True for column vector input.
%
%   ISCOLVECTOR(X) returns 1 if X is a column vector and 0 otherwise.
%
%   An array is considered a colum vector if the length along each
%   dimension, except possibly the first, is never larger than one.  That
%   is, SIZE(X, DIM) <= 1 for all DIM except possibly when DIM = 1.

%   Author:      Peter J. Acklam
%   Time-stamp:  2002-03-03 13:51:06 +0100
%   E-mail:      pjacklam@online.no
%   URL:         http://home.online.no/~pjacklam

   t = ndims(x) == 2 & size(x, 2) == 1;
