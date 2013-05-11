function U = unit(A, varargin)
%UNIT Normalize a matrix or vector.
%
%   Divides A by its norm. Computes the norm using the MATLAB norm function.
%   Passes optional second argument directly to norm, to indicate which kind
%   of norm to use.
%
%   For matrices...
%     UNIT(X) uses the 2-norm of X.
%     UNIT(X,2) is the same as UNIT(X).
%     UNIT(X,1) uses the 1-norm of X.
%     UNIT(X,inf) uses the infinity norm of X.
%     UNIT(X,'fro') uses the Frobenius norm of X.
%
%   For vectors...
%     UNIT(V,P) = V / sum(abs(V).^P)^(1/P).
%     UNIT(V) = V / norm(V,2).
%     UNIT(V,inf) = V / max(abs(V)).
%     UNIT(V,-inf) = V / min(abs(V)).
%
%   See also norm.

    U = A / norm(A, varargin{:});
end
