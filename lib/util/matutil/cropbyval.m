function y = cropbyval(x, val)
%CROPBYVAL Crop an array by
%
%   Y = CROPBYVAL(X, VAL) returns a subarray identical to X but where any
%   values VAL on the borders of X have been removed.
%
%   Y = CROPBYVAL(X) will determine VAL automatically by looking at the
%   borders of X and picking the most common value.  En error is given if
%   there is no value which is more common that the other values on the
%   border.  Note that automatically chosing VAL slowes down CROPBYVAL.
%
%   For example, if
%
%     X = [ 0  0  0  0
%           1  2  3  0
%           4  5  6  0 ]
%
%   then both Y = CROPBYVAL(X, 0) and Y = CROPBYVAL(X) will return
%
%     Y = [ 1  2  3
%           4  5  6 ]

%   Author:      Peter J. Acklam
%   Time-stamp:  2002-03-03 13:51:33 +0100
%   E-mail:      pjacklam@online.no
%   URL:         http://home.online.no/~pjacklam

   % check number of input arguments
   error(nargchk(1, 2, nargin));

   % if VAL is give, check it, otherwise compute default value
   if nargin == 2
      if any(size(val) ~= 1)
         error('VAL must be a scalar.');
      end
   else
      vals = [];
      sx = size(x);
      dx = ndims(x);
      c = cell(dx, 1);
      for i = 1:dx
         c{i} = 2:sx(i)-1;
      end
      t(prod(sx)) = logical(uint8(1));          % use uint8 to save memory
      t = reshape(t, sx);
      t(:) = t(end);
      t(c{:}) = 0;
      [lens, vals] = rlencode(sort(x(t(:))));
      len = max(lens);
      k = find(len == lens);
      if length(k) > 1
         error('No unique most common value on boundary.');
      end
      val = vals(k);
   end

   % get linear index values of elements different from VAL
   if isnan(val)
      i = find(~isnan(x));
   elseif val == 0
      i = find(x);
   else
      i = find(x ~= val);
   end

   % quick exit if output will be empty
   if isempty(i)
      y = [];
      return
   end

   % The following is based on code by
   % Doug Schwarz <douglas.schwarz@kodak.com>

   sx     = size(x);            % size vector of x
   dx     = length(sx);         % number of dimensions of x
   c      = cell(dx, 1);        % initialize list of subscript
   [c{:}] = ind2sub(sx, i);     % compute the subscripts
   for j = 1:dx
      c{j} = min(c{j}) : max(c{j});
   end
   y = x(c{:});
