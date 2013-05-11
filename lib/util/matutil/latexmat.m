function str = latexmat(mat, format)
%LATEXMAT Generate LaTeX code for a matrix.
%
%   STR = LATEXMAT(MAT, FORMAT) return the LaTeX code for the matrix MAT, with
%   the given FORMAT.
%
%   See HELP SPRINTF for more details about the FORMAT parameter.

%   Author:      Peter J. Acklam
%   Time-stamp:  2003-07-12 22:14:14 +0200
%   E-mail:      pjacklam@online.no
%   URL:         http://home.online.no/~pjacklam

   error(nargchk(2, 2, nargin));

   if ischar(mat)
      error('First argument can not be a string.');
   end

   if ~ischar(format)
      error('Second argument must be a string.');
   end

   [ r, c ] = size(mat);

   newline = sprintf('\n');

   str = [ '\left[ \begin{array}{' ...
               char(abs('c')*ones(1,c)) '}' newline ];

   for i = 1:r
      str = [ str ' ' ];
      for j = 1:c
         t = sprintf(format, real(mat(i,j)));
         if (imag(mat(i,j)) > 0)
            t = [ t '+' sprintf(format, imag(mat(i,j))) 'i' ];
         elseif (imag(mat(i,j)) < 0)
            t = [ t '-' sprintf(format, -imag(mat(i,j))) 'i' ];
         end
         str = [ str ' ' t ];
         if j < c
            str = [ str ' &' ];
         else
            if i < r
               str = [ str ' \\' newline ];
            else
               str = [ str newline ];
            end
         end
      end
   end
   str = [ str '\end{array} \right]' newline ];
