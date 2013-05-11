function b = isposint(x)
%ISPOSINT True for positive integers.
%
%   ISPOSINT(X) returns 1's where the elements of X are positive
%   integers and 0's where they are not.
%   For example, ISPOSINT([ 2e9 pi 3-5i -Inf ]) is [ 1 0 1 0 ].

%   Author:      Peter J. Acklam
%   Time-stamp:  2002-03-03 13:50:59 +0100
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
            b = x == round(x) & x > 0;
         case 'single'
            % convert input to double, compare, and convert back
            d = double(x);
            b = single(d == round(d)) & x > 0;
         case {'int8', 'int16', 'int32'}
            b = x > 0;
         case {'uint8', 'uint16', 'uint32'}
            % return an array of ones of the same class as input argument
            b = logical(repmat(feval(cls, 1), size(b)));
         otherwise
            error('Argument is of unrecognized class.');
         end
      end
   end
