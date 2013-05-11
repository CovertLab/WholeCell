function B = reparray(A, rep, siz)
%REPARRAY Replicate and tile an array or all subarrays of an array.
%
%   REPARRAY(A, REP) replicates the array A.  REP is a vector of integers
%   describing the number of replications to perform along each dimension.
%   Used this way, REPARRAY behaves like REPMAT.
%
%   REPARRAY(A, REP, SIZ) replicates the array A or all subarrays of A.  REP
%   is a vector of integers describing the number of replications to perform
%   along each dimension.  SIZ is a vector of integers describing the size
%   of each subarray.  An error will be given if A can not be divided into
%   subarrays of size SIZ.  If SIZ is omitted, SIZE(A) is used (see above).
%
%   In other words, REPARRAY(A, [M N ...], [P Q ...]) divides A into
%   P-by-Q-by-... subarrays and replicates each subarray M times along
%   the first dimension, N times along the second dimension, etc., so
%   the output is an M-by-N-by-... block array of P-by-Q-by-... arrays.
%
%   The size of the output array is only determined by the first two
%   arguments.  More precicely, REPARRAY(A, [M N ...], ...) returns an
%   array of size [M*SIZE(A,1) N*SIZE(A,2) ...].
%
%   Two special cases:
%
%      REPARRAY(A, [M N ...]) and REPARRAY(A, [M N ...], SIZE(A)) will
%      replicate and tile the whole array A.  Hence, it is equivalent to
%      REPMAT(A, [M N ...]).
%
%      REPARRAY(A, [M N ...], 1) will replicate each element of A.  When
%      A is a matrix of class double, REPARRAY(A, [M N], 1) is
%      equivalent to KRON(A,ONES(M,N)), but the former is generally a
%      lot faster since it does no multiplications.
%
%   Examples:
%
%       reparray(NaN, [2 3])
%       reparray(magic(2), [2 3], [1 1])
%       reparray(magic(2), [2 3], [1 2])
%       reparray(magic(2), [2 3], [2 1])
%       reparray(magic(2), [2 3], [2 2])    % = repmat(magic(2), [2 3])
%
%   See also REPMAT, MESHGRID, NDGRID.

%   Author:      Peter J. Acklam
%   Time-stamp:  2001-06-04 23:39:53 +0200
%   E-mail:      pjacklam@online.no
%   URL:         http://home.online.no/~pjacklam

   nargsin = nargin;
   nargchk(2, 3, nargin);
   if nargsin < 3
      siz = size(A);
   end

   if length(rep) == 1, rep = [rep rep]; end
   if length(siz) == 1, siz = [siz siz]; end

   if length(A) == 1 & all(siz == 1)

      nelems = prod(rep);
      if nelems > 0
         % Since B doesn't exist, the first statement creates a B with
         % the right size and type.  Then use scalar expansion to
         % fill the array.  Finally reshape to the specified size.
         B(nelems) = A;
         B(:) = A;
         B = reshape(B,rep);
      else
         B = A(ones(rep));
      end

   else

      Asiz = size(A);
      Adim = length(Asiz);
      Rdim = length(rep);
      Sdim = length(siz);
      Bdim = max([Adim Rdim Sdim]);

      Asiz = [Asiz ones(1,Bdim-Adim)];
      siz  = [siz  ones(1,Bdim-Sdim)];

      rat = Asiz./siz;
      if any(rat ~= round(rat))
         error('A can not be divided into subarrays as specified.');
      end

      subs{Bdim} = ':';
      subs(:) = subs(end);
      for i = 1:Rdim
         if rep(i) ~= 1
            ind = reshape(1:Asiz(i), [siz(i) 1 rat(i)]);
            subs{i} = ind(:,ones(rep(i),1),:);
         end
      end
      B = A(subs{:});

   end
