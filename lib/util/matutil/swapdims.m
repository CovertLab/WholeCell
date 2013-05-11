function y = swapdims(x, dim)
%SWAPDIMS Swap two dimensions.
%
%   SWAPDIMS(X, [DIM1 DIM2]) swaps dimensions DIM1 and DIM2.
%
%   SWAPDIMS(X, DIM) swaps dimensions 1 and DIM.
%
%   See also PERMUTE, IPERMUTE.

%   Author:      Peter J. Acklam
%   Time-stamp:  2003-04-17 10:37:11 +0200
%   E-mail:      pjacklam@online.no
%   URL:         http://home.online.no/~pjacklam

   error(nargchk(2, 2, nargin));

   dimlen = length(dim);
   if dimlen < 1 | dimlen > 2
      error('Second argument must be a vector with one or two elements.');
   else
      if any(dim < 1) | any(dim ~= round(dim))
         error('Dimensions must be positive integers.');
      end
      dim1 = dim(1);
      if dimlen < 2
         dim2 = 1;
      else
         dim2 = dim(2);
      end
   end

   if dim1 == dim2
      % both dimensions actually are the same dimension, so swapping
      % dimensions is a null-operation
      y = x;
      return
   end

   sdim1 = size(x, dim1);               % length of x along dimension 1
   sdim2 = size(x, dim2);               % length of x along dimension 2

   if sdim1 == 1 & sdim2 == 1
      % both dimensions are singleton, so swapping dimensions is a
      % null-operation
      y = x;
      return
   end

   nd = max([dim ndims(x)]);            % maximum dimension to deal with
   pv = 1:nd;                           % dimension permutation vector
   pv([dim1 dim2]) = [dim2 dim1];       % swap elements in permutation vector
   y = permute(x, pv);                  % now swap the dimensions in x
