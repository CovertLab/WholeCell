function varargout = gsexpand(varargin)
%GSEXPAND Generalized scalar expansion.
%
%   [XE, YE, ZE, ...] = GSEXPAND(X, Y, Z, ...) performs a generalized scalar
%   expansion on the input arguments.  The output arguments will all have
%   the same size.
%
%   All arguments are expanded (replicated) along singleton dimensions to
%   match the size of the other arguments.
%
%   See also RESIZE, SEXPAND.

%   Author:      Peter J. Acklam
%   Time-stamp:  2002-06-06 16:13:35 +0200
%   E-mail:      pjacklam@online.no
%   URL:         http://home.online.no/~pjacklam

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % check number of input and output arguments
   %

   nargsin  = nargin;
   nargsout = nargout;

   % check the number of input arguments
   if nargsin < 2
      error('Not enough input arguments.');
   end

   % check the number of output arguments
   if nargsout == 0
      % MATLAB convention: when a function is called with no output
      % arguments, return one output argument
      nargsout = 1;
   else
      if nargsout > nargsin
         error('Too many output arguments.');
      end
   end

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % find the size of the output arguments
   %

   % initialize common size vector and number of dimensions to that of a
   % scalar
   csize  = [1, 1];
   cdims  = 2;

   % iterate over the input arguments
   for i = 1:nargsin

      % size and number of dimensions for i'th input argument
      isize = size(varargin{i});
      idims = length(isize);

      if idims <= cdims
         % i'th argument has no more dimensions than the output will have;
         % only check the `idims' lowest dimensions
         dims = idims;
      else
         % i'th argument has more dimensions than the current common size
         % vector; only compare the `cdims' lowest dimensions
         dims = cdims;

         % update `cdims' for the next round
         cdims = idims;
      end

      % find non-singleton dimensions
      insdims = isize(1:dims) ~= 1;
      cnsdims = csize(1:dims) ~= 1;

      % check size compatibility
      if any( insdims & cnsdims & (isize(1:dims) ~= csize(1:dims)) )
         error('Lengths must match along non-singleton dimensions.');
      end

      csize(insdims) = isize(insdims);

   end

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % initialize output argument list and call other program
   %

   varargout = cell(1, nargsout);
   [varargout{:}] = resize(csize, varargin{:});
