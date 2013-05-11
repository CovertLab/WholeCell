function B = blockrot90(A, p, q, k)
%BLOCKROT90 Perform 90 degree block rotation on arrays.
%
%   BLOCKROT90(A,P,Q) is the 90 degree anticlockwise rotation of all P-by-Q
%   submatrices of the matrix A.
%
%   BLOCKROT90(A,P,Q,K) is the K*90 degree rotation of all P-by-Q submatrices
%   of A, K = +-1,+-2,...
%
%   BLOCKROT90 is vectorized along higher dimensions, so the same operations is
%   performed on all 2-D slices A(:,:,I,J,K,...).
%
%   See also ROT90.

%   Author:      Peter J. Acklam
%   Time-stamp:  2003-10-13 13:25:20 +0200
%   E-mail:      pjacklam@online.no
%   URL:         http://home.online.no/~pjacklam

   nargsin = nargin;
   error(nargchk(3, 4, nargsin));

   z = size(A);
   if nargsin < 4
      k = 1;
   else
      if any(size(k) ~= 1)
         error('K must be a scalar.');
      end
      k = k - 4*floor(k/4);     % map k to {0,1,2,3}
   end

   r = z(1)/p;
   s = z(2)/q;
   if (r ~= round(r)) | (s ~= round(s))
      error('A can not be subdivided into P-by-Q submatrices.');
   end

   u = [ p r q prod(z)/(z(1)*q) ];
   v = [ q*r p*s z(3:end) ];

   switch k
      case 1                    % 90 degrees anticlockwise
         B = reshape(A, u);
         B = reshape(permute(B(:,:,q:-1:1,:), [ 3 2 1 4 ]), v);
      case 2                    % 180 degrees
         B = reshape(A, u);
         B = reshape(B(p:-1:1,:,q:-1:1,:), z);
      case 3                    % 90 degrees clockwise
         B = reshape(A, u);
         B = reshape(permute(B(p:-1:1,:,:,:), [ 3 2 1 4 ]), v);
      otherwise                 % 0 degrees
         B = A;
   end
