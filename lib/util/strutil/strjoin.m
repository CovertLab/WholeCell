function str = strjoin(sep, varargin)
%STRJOIN Join strings in a cell array.
%
%   STRJOIN(SEP, STR1, STR2, ...) joins the separate strings STR1, STR2, ...
%   into a single string with fields separated by SEP, and returns that new
%   string.

%   Examples:
%
%     strjoin('-by-', '2', '3', '4')
%
%   returns '2-by-3-by-4'.
%
%     list = {'fee', 'fie', 'foe.m'};
%     strjoin('/', list{:}).
%
%   returns 'fee/fie/foe.m'.
%
%   This function is inspired by Perl' function join().

%   Author:      Peter J. Acklam
%   Time-stamp:  2003-10-13 11:13:55 +0200
%   E-mail:      pjacklam@online.no
%   URL:         http://home.online.no/~pjacklam

   % Check number of input arguments.
   error(nargchk(1, Inf, nargin));

   % Quick exit if output will be empty.
   if nargin == 1
      str = '';
      return
   end

   if isempty(sep)
      % special case: empty separator so use simple string concatenation
      str = [ varargin{:} ];
   else
      % varargin is a row vector, so fill second column with separator (using scalar
      % expansion) and concatenate but strip last separator
      varargin(2,:) = { sep };
      str = [ varargin{1:end-1} ];
   end
