function y = movsum(x, m)
%MOVSUM Moving sum.
%
%   Y = MOVSUM(X, M) returns a vector Y with the moving sum of length M, so
%   that
%
%     Y(1) = sum(X(1:M))
%     Y(2) = sum(X(2:M+1))
%     Y(3) = sum(X(3:M+2))
%     ...
%     Y(n-m+1) = sum(X(N-M+1:N))
%
%   To get the moving average, use MOVSUM(X, M)/M.

%   Author:      Peter J. Acklam
%   Time-stamp:  2003-08-17 13:33:09 +0200
%   E-mail:      pjacklam@online.no
%   URL:         http://home.online.no/~pjacklam

y = filter(ones(m,1), 1, x);
y = y(m:end);

return

sx = size(x);
k = find(sx > 1);
if length(k) > 1
   error('Input must be a vector.');
end
n = length(x);

i = 0 : n-m;
i = i(ones(m,1),:);
j = (1:m).';
j = j(:,ones(n-m+1,1));
y = sum(x(i+j), 1);

sy = sx;
sy(k) = n-m+1;
y = reshape(y, sy);
