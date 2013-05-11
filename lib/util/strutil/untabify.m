function newstr = untabify(str, n)
%UNTABIFY Convert tabs to spaces.
%
%   NEWSTR = UNTABIFY(STR) converts each tab character in STR to the
%   appropriate number of space characters.

%   Author:      Peter J. Acklam
%   Time-stamp:  2003-10-13 11:17:52 +0200
%   E-mail:      pjacklam@online.no
%   URL:         http://home.online.no/~pjacklam

   error(nargchk(1, 2, nargin));
   if ~ischar(str)
      error('First argument must be a string.');
   end
   if nargin < 2
      n = 8;
   else
      if any(size(n) ~= 1) | ~isnumeric(n) | n ~= round(n) | n < 0;
         error('Second argument must be a non-negative integer.');
      end
   end

   tab = sprintf('\t');                 % tab character

   strlen = length(str);                % length of input string

   newstr = '';
   if strlen == 0
      return
   end

   tabpos = find(str == tab);           % position of tabs in string
   tabpos = [ tabpos strlen+1 ];        % add "end of line mark"

   ntabs = length(tabpos);              % number of tabs

   newstr = str(1:tabpos(1)-1);
   newstrlen = length(newstr);

   for i = 1:ntabs-1
      nspc = n*ceil((newstrlen+1)/n) - newstrlen;
      newstr = [ newstr blanks(nspc) str(tabpos(i)+1:tabpos(i+1)-1) ];
      newstrlen = length(newstr);
   end
