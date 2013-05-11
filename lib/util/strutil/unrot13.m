function varargout = unrot13(varargin)
%UNROT13 Rot13 (Caesar) decryption.
%
%   UNROT13(STR1, STR2, ...) applies a rot13 (Caesar) decryption on all
%   input strings and returns the same number of strings.
%
%   For instance, unrot13('Uryyb Jbeyq!') returns 'Hello World!'.

%   Author:      Peter J. Acklam
%   Time-stamp:  2003-10-13 11:17:40 +0200
%   E-mail:      pjacklam@online.no
%   URL:         http://home.online.no/~pjacklam

   % check number of input arguments
   error(nargchk(1, Inf, nargin));

   % use the fact that a double rot13 operation is a no-op
   varargout = cell(size(varargin));
   [varargout{:}] = rot13(varargin{:});
