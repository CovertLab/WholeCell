function dim = fnsdim(x)
%FNSDIM First non-singleton dimension.
%
%   FNDSDIM(X) returns the first non-singleton dimension of X.  If X has no
%   non-singleton dimensions, then 1 is returned.
%
%   A singleton dimension is a dimension along which the length is one.
%
%   There are numerous MATLAB functions that default to operate along the
%   first non-singleton dimension.  This function is useful for finding this
%   default dimension.
%
%   See also SIZE, NDIMS.

%   Author:      Peter J. Acklam
%   Time-stamp:  2001-07-07 22:56:52 +0200
%   E-mail:      pjacklam@online.no
%   URL:         http://home.online.no/~pjacklam

   error(nargchk(1, 1, nargin));

   k = find(size(x) ~= 1);
   if isempty(k)
      dim = 1;
   else
      dim = k(1);
   end
