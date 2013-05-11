function varargout = ndimgrid(varargin)
%NDIMGRID Generation of arrays for N-D functions and interpolation.
%
%   [X1,X2,...,XN] = NDIMGRID(x1,x2,...,xN) transforms the domain specified
%   by vectors x1,x2,...,xN into arrays X1,X2,...,XN that can be used for
%   the evaluation of functions of N variables and N-D interpolation.  The
%   i-th dimension of the output array Xi are copies of elements of the
%   vector xi.
%
%   Differences between NDIMGRID and NDGRID:
%
%   * NDIMGRID works on any array.  NDGRID only works on numerical arrays.
%
%   * NDIMGRID expands the input arguments strictly according to the number
%     of input arguments.  As a consequence, if input is a single vector,
%     NDIMGRID returns a vector, whereas NDGRID expands it to a matrix as if
%     two input arguments were given.
%
%   * NDIMGRID gives an error if too many output arguments are given,
%     whereas NDGRID returns the extra output arguments as empty arrays.
%
%   * NDIMGRID is faster than NDGRID.
%
%   See also NDGRID, MESHGRID, INTERPN.

   % check number of input arguments
   nin = nargin;
   error(nargchk(1, Inf, nin));

   nout = max(nargout, nin);

   for i = nin : -1 : 1
      numel(i) = prod(size(varargin{i}));   % num of elems in i'th arg
      subs{i}  = ones(numel(i),1);          % add subscript cell array entry
   end

   varargout = cell(1, nout);           % initialize output arg list
   one = ones(1, max(nout, 2));         % initialize reshape size vector

   for i = 1 : min(nout,nin)
      % make sure numerical arrays are full
      if isnumeric(varargin{i})
         varargin{i} = full(varargin{i});
      end
      % reshape i'th argument into a vector x
      s = one;                          % size(x, j) = 1 for all j ~= i
      s(i) = numel(i);                  % size(x, i) = numel(i)
      x = reshape(varargin{i}, s);      % create x
      % replicate the vector to create an array
      s = subs;                         % replicate along all dimensions...
      s{i} = ':';                       % ...except the i'th dimension
      varargout{i} = x(s{:});
   end
