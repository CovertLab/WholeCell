function b = iseven(x)
%ISEVEN True for even numbers.
%
%   ISEVEN(X) returns 1's where the elements of X are even numbers and 0's
%   where they are not.
%
%   See also ISINT, ISODD.

%   Author:      Peter J. Acklam
%   Time-stamp:  2003-11-03 18:55:21 +0100
%   E-mail:      pjacklam@online.no
%   URL:         http://home.online.no/~pjacklam

   error(nargchk(1, 1, nargin));
   if ~isnumeric(x)
      error('Argument must be a numeric array.');
   end

   cls = class(x);                      % class of input argument
   if isempty(x)
      b = feval(cls, x);                % return empty array of same class
   else
      switch cls
         case 'double'
            b = ~mod(x, 2);
         case 'single'
            % "mod" is not defined for class "single"; so convert input to double, compare,
            % and convert back
            b = single(~mod(double(x), 2));
         case {'uint8', 'uint16', 'uint32', 'uint64'}
            b = ~bitand(x, 1);
         case {'int8', 'int16', 'int32', 'int64'}
            error('Not implemented for classes int8, int16, int32, and int64.');
         otherwise
            error('Argument is of unrecognized class.');
      end
   end
