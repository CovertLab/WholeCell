function out = cplexoptimget (options, name, default)
%%
% Purpose
% Retrieve optimization options values.
%
% Syntax
%    out = cplexoptimget(options,'param')
%    out = cplexoptimget(options,'param',default)
%
% Description
% val = cplexoptimget(options,'param') returns the value of the specified
% option in the optimization options structure options.
%
% val = cplexoptimget(options,'param',default) returns default if the
% specified option is not defined in the optimization options structure
% options.
%
%  See also cplexoptimset
%

% ---------------------------------------------------------------------------
% File: cplexoptimget.m
% ---------------------------------------------------------------------------
% Licensed Materials - Property of IBM
% 5725-A06 5725-A29 5724-Y48 5724-Y49 5724-Y54 5724-Y55
% Copyright IBM Corporation 2008, 2010. All Rights Reserved.
%
% US Government Users Restricted Rights - Use, duplication or
% disclosure restricted by GSA ADP Schedule Contract with
% IBM Corp.
% ---------------------------------------------------------------------------
