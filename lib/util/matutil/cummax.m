function y = cummax(x, dim)
%CUMMAX Cumulative maximum of elements.
%
%   For vectors, CUMMAX(X) is a vector containing the cumulative maximum of
%   the elements of X, that is X(i) is the MAX(X(1:i)).
%
%   For matrices, CUMMAX(X) is a matrix the same size as X containing the
%   cumulative maximum over each column.
%
%   For N-D arrays, CUMMAX(X) operates along the first non-singleton
%   dimension.
%
%   CUMMAX(X, DIM) works along the dimension DIM.
%
%   For small vectors, it is much faster to use the following solution by
%   Ruud Brekelmans <R.C.M.Brekelmans@kub.nl>:
%
%     y = max(x(cumsum(triu(ones(length(x))))));
%
%   See also CUMSUM, CUMPROD.

%   Author:      Peter J. Acklam
%   Time-stamp:  2001-06-02 00:30:03 +0200
%   E-mail:      pjacklam@online.no
%   URL:         http://home.online.no/~pjacklam

   nargsin = nargin;
   error(nargchk(1, 2, nargsin));

   s = size(x);
   if nargsin < 2               % if DIM is not specified
      dim = find(s ~= 1);       % find non-singleton dimensions
      if isempty(dim)           % if only singleton dimensions
         y = x;                 %   x is scalar so output = input
         return;                %   bail out
      else                      % at least one non-singleton dimension
         dim = dim(1);          %   so pick first non-singleton dimension
      end
   else                         % DIM is specified
      if any(size(dim) ~= 1) | ~isreal(dim) | dim ~= round(dim) | dim < 1
         error('DIM must be a real scalar positive integer.');
      end
   end

   % output = input if length along dimension <= 1
   if size(x, dim) <= 1
      y = x;
      return;
   end

   % it is easiest to work with the first dimension of a 2D array
   v = 1:ndims(x);              % initialize dimension permutation vector
   v([1 dim]) = [dim 1];        % we want to swap dimensions 1 and DIM
   x = permute(x, v);           % now do the swapping
   x = x(:,:);                  % convert from ND to 2D

   m = s(v(1)) - 1;

   % here is how this works:  1) find each value that is larger than the
   % next value 2) if none was found, bail out 3) adjust the linear indices
   % because of the fact that the matrices `x(1:end-1,:)' and `x(2:end,:)'
   % has one fewer row than `x' 4) in all the cases where a value is larger
   % than the next, let the next be identical to the previous
   while 1
      i = find(x(1:end-1,:) > x(2:end,:));  % find values > next values
      if isempty(i), break; end;            % bail out if we're done
      i = i + floor((i-1)/m);               % adjust linear indices
      x(i+1) = x(i);                        % let next = previous
   end

   y = reshape(x, s(v));        % convert from 2D to ND
   y = permute(y, v);           % swap dimensions 1 and DIM again
