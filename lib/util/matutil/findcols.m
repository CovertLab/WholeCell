function k = findcols(A, b)
%FINDCOLS Find indices of a given column within a matrix.
%
%   FINDCOLS(A, B) returns a row vector with the indices of the columns
%   in the matrix A that are identical to the column vector B.  If no
%   columns in A are identical to B, an empty vector is returned.
%
%   The methods uses a for-loop, but it uses less memory and is in many
%   cases a lot faster than the vectorized methods
%
%      find( all( A == repmat(b, 1, size(A, 2)), 1 ) )
%      find( all( A == b(:,ones(size(A, 2), 1)), 1 ) )
%
%   See also FIND, FINDROWS.

%   Author:      Peter J. Acklam
%   Time-stamp:  2002-03-03 13:51:19 +0100
%   E-mail:      pjacklam@online.no
%   URL:         http://home.online.no/~pjacklam

   k = find( A(1,:) == b(1) );
   for j = 2:size(A, 1)
      k = k( A(j,k) == b(j) );
      if isempty(k)
         return
      end
   end
