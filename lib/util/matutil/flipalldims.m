function y = flipalldims(x)
%FLIPALLDIMS Flip array along all dimensions.
%
%   FLIPALLDIMS(X) flips the array X along all dimension.
%
%   See also FLIPDIM.

%   Author:      Peter J. Acklam
%   Time-stamp:  2001-06-04 23:42:03 +0200
%   E-mail:      pjacklam@online.no
%   URL:         http://home.online.no/~pjacklam

   % check number of input arguments
   error(nargchk(1, 1, nargin));

   % quick exit if input is empty
   if isempty(x)
      y = x;
      return;
   end

   % get size of X and find non-singleton dimensions
   sx   = size(x);
   dims = find(sx > 1);

   % initialize subscript cell array
   idx{ndims(x)} = ':';
   idx(:) = idx(end);

   % flip along all non-singleton dimensions (flipping along singleton
   % dimensions is a null-operation)
   for i = 1:length(dims)
      dim = dims(i);
      idx{dim} = sx(dim):-1:1;
   end

   y = x(idx{:});
