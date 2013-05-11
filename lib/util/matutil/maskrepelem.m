function y = maskrepelem(x, mask)
%MASKREPELEM Replicate elements according to a mask array.
%
%   MASKREPELEM(X, MASK) replicates each element of X so it becomes a block
%   with the size of MASK and the element inserted at the positions where
%   MASK has a non-zero value.  X must be a matrix, but MASK may be any ND
%   array.  X and MASK may be of any class.
%
%   When X and MASK are double arrays, and MASK is a matrix of zeros and
%   ones, MASKREPELEM(X, MASK) returns the same as KRON(X, MASK).
%
%   For example, maskrepelem(magic(2), [0 1;1 0]) returns
%
%      [ 0  1  0  3
%        1  0  3  0
%        0  4  0  2
%        4  0  2  0 ]

%   Author:      Peter J. Acklam
%   Time-stamp:  2002-03-03 13:20:49 +0100
%   E-mail:      pjacklam@online.no
%   URL:         http://home.online.no/~pjacklam

% some error checking
error(nargchk(2, 2, nargin));
if ndims(x) > 2
   error('X must be a 2D matrix.');
end

% get sizes, number of elements and number of dimensions
sx = size(x);                       % size of x
nx = prod(sx);                      % number of elements in x
sm = size(mask);                    % size of mask
nm = prod(sm);                      % number of elements in mask
dm = ndims(mask);                   % number of dimensions in mask

% initialize output; replicate and insert the elements of x
i = find(mask);                     % index of non-zero elements
x = reshape(x, 1, nx);              % make x a row vector
y(nm, nx) = feval(class(x), 0);     % let y have same class as x
y(i,:) = x(ones(length(i),1),:);    % replicate and insert elements

% manipulate the blocks
y = reshape(y, [ sm sx ]);
y = permute(y, [ 1 1+dm 2 2+dm 3:dm ]);
y = reshape(y, [ sx(1)*sm(1) sx(2)*sm(2) sm(3:end) ]);
