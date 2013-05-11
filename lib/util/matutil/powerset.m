function X = powerset(varargin)
%POWERSET All combinations of elements taken from different vectors.
%
%   X = POWERSET(A, B, C, ...) returns a matrix containing all combinations
%   of elements taken one at a time from the input vectors.
%
%     [ A(1)   B(1)    C(1)  ]
%     [ A(1)   B(1)    C(2)  ]
%              ...
%     [A(end) B(end) C(end-1)]
%     [A(end) B(end)  C(end) ]
%
%   The number or rows in X is LENGTH(A)*LENGTH(B)*... and the number of
%   columns is the same as the number of input argument.

%   Author:      Peter J. Acklam
%   Time-stamp:  2002-03-03 13:52:36 +0100
%   E-mail:      pjacklam@online.no
%   URL:         http://home.online.no/~pjacklam

   n = length(varargin);
   tmp = cell(1, n);
   [tmp{n:-1:1}] = ndimgrid(varargin{n:-1:1});
   X = reshape(cat(n+1, tmp{:}), [prod(size(tmp{1})) n]);
