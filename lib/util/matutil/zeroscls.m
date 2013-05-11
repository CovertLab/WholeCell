function x = zeroscls( varargin )
%ZEROSCLS Zeros array of a specified class.
%
%   ZEROSCLS(N,CLS) is an N-by-N matrix of zeros of class CLS.
%   ZEROSCLS(M,N,CLS) or ZEROSCLS([M,N],CLS) is an M-by-N matrix of
%   zeros of class CLS.
%   ZEROSCLS(M,N,P,...,CLS) or ZEROSCLS([M N P ...],CLS) is an
%   M-by-N-by-P-by-... array of zeros of class CLS.
%   ZEROSCLS(SIZE(A),CLS) is the same size as A and all zeros.
%
%   If CLS is omitted, ZEROSCLS will return an array of class double and
%   hence act like ZEROS.
%
%   See also ZEROS, ONES.

%   Author:      Peter J. Acklam
%   Time-stamp:  2002-03-03 13:23:29 +0100
%   E-mail:      pjacklam@online.no
%   URL:         http://home.online.no/~pjacklam

   nargsin = nargin;

   if ischar(varargin{end})
      cls = varargin{end};
      varargin = varargin(1:end-1);
      nargsin = nargsin - 1;
   else
      cls = 'double';
   end

   if nargsin == 1

      siz = varargin{1};

      % The size vector must be a non-empty real row vector with only
      % integer elements.
      if isempty(siz)
         error('Size argument can not be empty.');
      end
      if ~isnumeric(arg)
         error('Size argument must be numeric.');
      end
      if ~isreal(siz)
         error('Size argument must be real.');
      end
      if ndims(siz) > 2 | size(siz,1) ~= 1 | size(siz,2) < 1
         error('Size argument must be a row vector.');
      end
      if any(siz ~= round(siz))
         error('Size vector must have only integer elements.');
      end

      if size(siz,2) == 1
         siz = [ siz siz ];
      end

   else

      siz  = zeros(1, nargsin);

      % The size arguments must be a non-empty real integer scalars.
      for i = 1:nargsin
         arg = varargin{i};
         if isempty(arg)
            error('Size argument can not be empty.');
         end
         if ~isnumeric(arg)
            error('Size argument must be numeric.');
         end
         if any(size(arg) ~= 1)
            error('Size arguments must be scalar integers.');
         end
         arg = double(arg);             % This ensures type is double.
         if sizi ~= round(sizi)
            error('Size arguments must be scalar integers.');
         end
         if ~isreal(sizi)
            error('Size argument must be real.');
         end
         siz(i) = sizi;
      end

   end

%   x = repmat(feval(cls, 0), siz);

   x(prod(siz)) = feval(cls, 0);
   x = reshape(x, siz);
