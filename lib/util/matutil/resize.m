function varargout = resize(csize, varargin)
%RESIZE Resize arguments to specified size.
%
%   [XE, YE, ZE, ...] = RESIZE(S, X, Y, Z, ...) resizes each of the input
%   arguments to the size S.  The output arguments will all have the same
%   size.  The resizing is done by truncation or recycling, as
%   necessary, along each dimension.  For instance,
%
%      resize([1 2], [7 8 9])  returns  [7 8]           % truncation
%      resize([1 3], [7 8 9])  returns  [7 8 9]         % no change
%      resize([1 5], [7 8 9])  returns  [7 8 9 7 8]     % recycling
%      resize([1 5], [7 8])    returns  [7 8 7 8 7]     % recycling
%      resize([1 5], 7)        returns  [7 7 7 7 7]     % recycling
%
%   When the size along the dimension to be replicated is one (the last case
%   above), the recycling degenerates to a duplication of elements.
%
%   When expanding an array X to the size S, the following applies to each
%   dimension I:
%
%     If SIZE(X, I) = S(I), then no change is made.
%
%     If SIZE(X, I) > S(I), then X is truncated along the dimension I to the
%     length S(I).
%
%     If SIZE(X, I) < S(I), then X will be expanded along dimension I by
%     recycling, i.e., X will be replicated along dimension I according to
%     an index vector [ 1 2 ... SIZE(X, I) 1 2 ... ] which has length S(I).
%
%   The special case when SIZE(X, I) == 1 and S(I) > 1 is called `scalar
%   expansion', i.e., X is simply duplicated S(I) times along dimension I.
%
%   See also SIZE, REPMAT.

%   Author:      Peter J. Acklam
%   Time-stamp:  2001-05-14 02:31:18 +0200
%   E-mail:      pjacklam@online.no
%   URL:         http://home.online.no/~pjacklam

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Check number of input arguments.
   %

   nargsin = nargin;
   if nargsin < 2
      error('Not enough input arguments.');
   end

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Check the size vector.  Apply the same requirements to the size vector
   % for RESIZE as the size vector for RESHAPE.
   %

   % Make sure it has at least two elements.
   vecsiz = size(csize);
   if prod(vecsiz) < 2
      error('Size vector must have at least two elements.');
   end

   % Make sure it is a row vector.  There must be only one non-singleton
   % dimension and that dimension must be the second.
   k = find(vecsiz ~= 1);               % find non-singleton dimensions
   if length(k) ~= 1 | k ~= 2
      error('Size vector must be a row vector with integer elements.');
   end

   % Make sure it contains only real numbers.
   if ~isnumeric(csize) | ~isreal(csize)
      error('Size vector must contain only real non-negative integers.');
   end

   % Make sure the numbers are positive integers.
   csize = double(csize);
   if any(isnan(csize)) | any(csize < 0) | any(csize ~= round(csize))
      error('Size vector must contain only real non-negative integers.');
   end

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Get number of dimensions in output arguments.
   %

   k = find(csize ~= 1);        % find non-singleton dimensions
   if isempty(k)
      cdims = 2;                % all arrays have at least two dimensions
   else
      cdims = k(end);           % last non-singleton dimension
      cdims = max(cdims, 2);    % must have at least two dimensions
   end
   csize = csize(1:cdims);      % remove any trailing singleton dimensions

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Prepare output argument list.  The number of arrays to resize should
   % depend on `nargout', not `nargin', to avoid resizing more arrays than
   % those that will be returned
   %

   nargsout = nargout;

   % When called with no output arguments, return one output argument.
   if nargsout == 0
      nargsout = 1;
   end

   if nargsout > length(varargin)
      error('More output arguments than input arguments specified.');
   end

   % Initialize output argument list.
   varargout = cell(1, nargsout);

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Define other variables that are used in the loop.
   %

   % number of elements in output array
   nelems = prod(csize);

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Special case: all output arguments are empty.  Still, we must reshape
   % the input argument to the specified size, rather than simply returning
   % `[]', because the input might not be of class `double' and it might not
   % be `0-by-0'.
   %

   if nelems == 0
      for i = 1 : nargsout
         varargout{i} = reshape(varargin{i}([]), csize);
      end
      return
   end

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % General case: resize all input arguments to some non-empty array.
   %

   % All input arrays must be non-empty if the output arrays are to be
   % non-empty.
   if any(cellfun('isempty', varargin))
      error('Can''t resize an empty array to a non-empty array.');
   end

   % Initialize default subscript array.
   subs{5} = ':';
   subs(:) = subs(end);

   for i = 1 : nargsout

      % Size of i'th input argument.
      isize = size(varargin{i});

      % If i'th input argument already has the correct size, just assign
      % directly.
      if isequal(csize, isize)
         varargout{i} = varargin{i};

      % If i'th input argument is a scalar; initialize, fill and reshape.
      elseif all(isize == 1)
         varargout{i}(nelems) = varargin{i};                    % initialize
         varargout{i}(:)      = varargin{i};                    % fill
         varargout{i}         = reshape(varargout{i}, csize);   % reshape

      % The general case: recycle or truncate as necessary.
      else
         % Number of dimensions of i'th input argument.
         idims = length(isize);

         % Initialize list of subscripts.
         subsi = subs;

         % Find the dimensions for which the lengths are different.
         k = find(isize(1:idims) ~= csize(1:idims));

         % Go through all dimensions for which the lengths are different.
         for j = 1 : length(k)

            % If length along j'th dimension is too long - truncate.
            if isize(k(j)) > csize(k(j))
               subsi{k(j)} = 1 : csize(k(j));

            % If length along j'th dimension is too short - recycle.
            elseif isize(k(j)) < csize(k(j))
               subsi{k(j)} = rem(0 : csize(k(j)) - 1, isize(k(j))) + 1;

            % We should never get here.
            else
               error('Internal error');
            end

         end

         % Now resize the i'th input argument.
         varargout{i} = varargin{i}(subsi{:});

      end
   end
