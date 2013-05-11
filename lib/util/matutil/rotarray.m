function B = rotarray(A, n, dim)
%ROTARRAY Rotate array or subarrays.
%
%   ROTARRAY(A) is the 90 degree anticlockwise rotation of matrix A.
%   ROTARRAY(A, N) is the K*90 degree rotation of A, N = 0,+-1,+-2,...
%
%   Example,
%      A = [ 1 2 3      B = rotarray(A) = [ 3 6
%            4 5 6 ]                        2 5
%                                           1 4 ]
%
%   If A is an ND array, the same operation is performed on all 2-D slices
%   A(:,:,I,J,K,...).
%
%   ROTARRAY(A, K, DIM), where DIM is a vector with two integers, will perform
%   the rotation anticlockwise around an axis perpendicular to a plane through
%   the specified dimensions.  When DIM is omitted, the dimension vector [1 2]
%   is used, which gives the same behaviour as ROT90.
%
%   ROTARRAY(A, K, DIM) where DIM is a scalar, is equivalent to
%   ROTARRAY(A, K, [DIM DIM+1]).
%
%   ROTARRAY(A, K, [DIM1 DIM2]) is equivalent to
%   ROTARRAY(A, 4-K, [DIM2 DIM1]).
%
%   See also ROT90.

%   Author:      Peter J. Acklam
%   Time-stamp:  2003-10-13 13:25:03 +0200
%   E-mail:      pjacklam@online.no
%   URL:         http://home.online.no/~pjacklam

   % check number of input arguments
   nargsin = nargin;
   error(nargchk(1, 3, nargsin));

   % get N
   if (nargsin < 2) | isempty(n)
      n = 1;
   else
      if any(size(n) ~= 1) | (n ~= round(n))
         error('N must be a scalar integer.');
      end
      n = n - 4*floor(n/4);     % map n to {0,1,2,3}
   end

   % get dimension vector DIM
   if nargsin < 3
      dim = [1 2];
   else
      nd = prod(size(dim));
      if nd < 1 | nd > 2 | any(dim <= 0) | any(dim ~= round(dim))
         error('DIM must contain one or two positive integers.');
      end
      if nd == 1
         dim = [dim dim+1];
      end
   end

   dim1  = dim(1);
   dim2  = dim(2);
   sdim1 = size(A, dim1);
   sdim2 = size(A, dim2);

   if (n == 0) | ((sdim1 <= 1) & (sdim2 <= 1))
      % special case when the rotation is a null-operation (output will be
      % identical to input)
      B = A;
   else
      % largest dimension number we have to deal with
      nd = max([ndims(A) dim1 dim2]);

      % initialize subscript cell array
      v{nd} = ':';
      v(:) = v(end);

      switch n
         case 1                 % 90 degrees anticlockwise
            v{dim2} = sdim2 : -1 : 1;
            d = 1:nd;
            d([ dim1 dim2 ]) = [ dim2 dim1 ];
            B = permute(A(v{:}), d);
         case 2                 % 180 degrees anticlockwise
            v{dim1} = sdim1 : -1 : 1;
            v{dim2} = sdim2 : -1 : 1;
            B = A(v{:});
         case 3                 % 270 degrees anticlockwise
            v{dim1} = sdim1 : -1 : 1;
            d = 1:nd;
            d([ dim1 dim2 ]) = [ dim2 dim1 ];
            B = permute(A(v{:}), d);
      end
   end
