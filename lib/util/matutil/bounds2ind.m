function idx = bounds2ind(lo, hi)
%BOUNDS2IND Create index vector given bounds for index values.
%
%   BOUNDS2IND(LO, HI) will, when given two vectors LO and HI, return the
%   index vector [LO(1):HI(1) LO(2):HI(2) ...].

%   Author:      Peter J. Acklam
%   Time-stamp:  2001-10-23 17:05:58 +0200
%   E-mail:      pjacklam@online.no
%   URL:         http://home.online.no/~pjacklam

   % check number of input arguments
   error(nargchk(2, 2, nargin));

   % keep only runs of positive length.
   i = lo <= hi;
   lo = lo(i);
   hi = hi(i);

   m   = length(lo);            % length of input vectors
   len = hi - lo + 1;           % length of each run
   n   = sum(len);              % length of output vector
   idx = ones(1, n);            % initialize output vector

   idx(1) = lo(1);
   len(1) = len(1) + 1;
   idx(cumsum(len(1:end-1))) = lo(2:m) - hi(1:m-1);
   idx = cumsum(idx);
