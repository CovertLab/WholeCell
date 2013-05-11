function B = mergearray(A)
%MERGEARRAY Merge all arrays in a cell array into a single array.
%
%   MERGEARRAY(A) merges all arrays in the cell array A into a single
%   array.  All arrays in A must have the same size.
%
%   See also CAT, SPLITARRAY.

%   Author:      Peter J. Acklam
%   Time-stamp:  2002-03-03 13:50:41 +0100
%   E-mail:      pjacklam@online.no
%   URL:         http://home.online.no/~pjacklam

error(nargchk(1, 1, nargin));

Asiz = size(A);                     % size of A
Adim = length(Asiz);                % dimensions in A
Ssiz = size(A{1});                  % size in each subarray
Sdim = length(Ssiz);                % dimensions in each subarray
Bdim = max(Adim, Sdim);             % dimensions in B (output array)
Asiz = [Asiz ones(1,Bdim-Adim)];    % size of A (padded)
Ssiz = [Ssiz ones(1,Bdim-Sdim)];    % size of each subarray (padded)
Bsiz = Ssiz .* Asiz;                % size of B (output array)

% B is now [ Ssiz(1) Ssiz(2) ... Ssiz(end) Asiz(1) Asiz(2) ... Ssiz(end)].
B = reshape(cat(Bdim+1, A{:}), [Ssiz Asiz]);

% A is now [ Ssiz(1) Asiz(1) Ssiz(2) Asiz(2) ... Ssiz(end) Asiz(end) ].
B = permute(B, reshape([1:Bdim ; Bdim+1:2*Bdim], [1 2*Bdim]));

% A is now [ Ssiz(1)*Asiz(1) Ssiz(2)*Asiz(2) ... Ssiz(end)*Asiz(end) ].
B = reshape(B, Bsiz);
