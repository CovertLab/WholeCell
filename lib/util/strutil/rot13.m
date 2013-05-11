function varargout = rot13(varargin)
%ROT13  Rot13 (Caesar) encryption.
%
%   ROT13(STR1, STR2, ...) applies a rot13 (Caesar) encryption on all input
%   strings and returns the same number of strings.
%
%   For instance, rot13('Hello World!') returns 'Uryyb Jbeyq!'.

%   Author:      Peter J. Acklam
%   Time-stamp:  2003-10-13 11:13:10 +0200
%   E-mail:      pjacklam@online.no
%   URL:         http://home.online.no/~pjacklam

   % check number of input arguments
   error(nargchk(1, Inf, nargin));

   % create mapping vectors
   i = [ 13:25, 0:12 ];
   j = [ 0:'A'-1,       ...     % untranformed
         'A'+i,         ...     % upper case letters
         'Z'+1:'a'-1,   ...     % untransformed
         'a'+i,         ...     % lower case letters
         'z'+1:255,     ...     % untransformed
       ];

   varargout = cell(size(varargin));
   for i = 1:nargin
      varargout{i} = char(j(1+varargin{i}));
   end
