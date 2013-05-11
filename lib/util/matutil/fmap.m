function array = fmap(fun, array, varargin)
%FMAP   Evaluate a function for each element of an array.
%
%   ARRAY2 = FMAP(FUN, ARRAY) evaluates the function FUN on all
%   elements of ARRAY.  FUN must be a string and ARRAY1 must be numeric
%   array, a character array or a cell array.
%
%   ARRAY2 = FMAP(FUN, ARRAY, X1, X2, ...) passes extra arguments to
%   FUN.
%
%   Examples:
%
%     fmap('sqrt', [ 4 9 ])               returns  [ 2 3 ]
%
%     fmap('sqrt', { 4 9 })               returns  { 2 3 }
%
%     fmap('power', { 4 9 }, 3)           returns  { 64 729 }
%
%     fmap('all', { [ 1 1 0 ] [ 1 1 ] })  returns  { 0 1 }

%   Author:      Peter J. Acklam
%   Time-stamp:  2002-03-03 13:51:14 +0100
%   E-mail:      pjacklam@online.no
%   URL:         http://home.online.no/~pjacklam

   error(nargchk(2, Inf, nargin));

   if isnumeric(array) | ischar(array)

%      eval(sprintf([ ...
%         'for i = 1 : prod(size(array))\n' ...
%         '   array(i) = %s(array(i), varargin{:});\n' ...
%         'end\n' ], fun));

      for i = 1 : prod(size(array))
         array(i) = feval(fun, array(i), varargin{:});
      end

   elseif isa(array, 'cell')

%      eval(sprintf([ ...
%         'for i = 1 : prod(size(array))\n' ...
%         '   array{i} = %s(array{i}, varargin{:});\n' ...
%         'end\n' ], fun));

      for i = 1 : prod(size(array))
         array{i} = feval(fun, array{i}, varargin{:});
      end

   else

      error(['No support for arrays of class "' class(array) '".' ]);

   end
