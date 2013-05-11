function d = distmat(x, y)
%DISTMAT Distance matrix.
%
%   D = DISTMAT(X, Y) returns the distance matrix with all distances
%   between the elements in X and Y. Distance is Euclidean.
%
%   DISTMAT(X) is the same as DISTMAT(X,X).
%
%   For matrix arguments, each row is considered to be a point.

%   Author:      Peter J. Acklam
%   Time-stamp:  2002-03-03 13:51:31 +0100
%   E-mail:      pjacklam@online.no
%   URL:         http://home.online.no/~pjacklam

   if (nargin == 2)
      d = distmat1(x, y);
   else
      d = distmat1(x);
   end
