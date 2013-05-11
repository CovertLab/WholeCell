function tf = strlexcmp(a, b)
%STRLEXCMP Lexicographic comparison of two strings.
%
%   STRLEXCMP(A, B) returns -1, 0, or 1 depending on whether the left argument
%   is stringwise less than, equal to, or greater than the right argument.
%
%   This is a MATLAB version of the Perl `cmp' operator.
%
%   See also EQ, ISEQUAL.

%   Author:      Peter J. Acklam
%   Time-stamp:  2003-10-13 11:16:31 +0200
%   E-mail:      pjacklam@online.no
%   URL:         http://home.online.no/~pjacklam

   % check arguments
   error(nargchk(2, 2, nargin));
   if ~ischar(a) | ~ischar(b)
      error('Both arguments must be char arrays (strings).');
   end

   % get lengths of strings
   na = length(a);
   nb = length(b);
   n = min(na, nb);

   % find characters that differ
   k = find(a(1:n) ~= b(1:n));
   if isempty(k)
      % all characters are identical -- compare lengths
      tf = sign(na - nb);
   else
      % compare first character that is different
      k = k(1);
      tf = sign(a(k) - b(k));
   end
