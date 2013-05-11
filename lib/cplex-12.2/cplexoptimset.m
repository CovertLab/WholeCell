function options = cplexoptimset (varargin)
%%
% Purpose
% Create or edit optimization options structure.
%
% Syntax
%    options = cplexoptimset('param1',value1,'param2',value2,...)
%    cplexoptimset
%    options = cplexoptimset
%    options = cplexoptimset('cplex')
%    options = cplexoptimset(oldopts,'param1',value1,...)
%    options = cplexoptimset(oldopts,newopts)
%
% Description
% options = cplexoptimset('param1',value1,'param2',value2,...) creates an
% optimization options structure called options, in which the specified
% options (param) have specified values. The specified options correspond
% to the options in the MATLAB(R) Optimization Toolbox.Any unspecified 
% options are set to [] (options with value [] indicate that the default 
% value for that option should be used when options is passed to the 
% optimization function).
%
% cplexoptimset with no input or output arguments displays a complete list
% of options with their valid values.
%
% options = cplexoptimset (with no input arguments) creates an option
% structure options where all fields are set to [], which instructs
% CPLEX(R) to use all default parameter values. 
%
% options = cplexoptimset ('cplex') creates a structure options, which 
% contains all of the CPLEX parameters.  Use this method when you need to 
% set parameters that do not correspond to the options in the MATLAB 
% Optimization Toolbox.
%
% options = cplexoptimset(oldopts,'param1',value1,...) creates a copy of
% oldopts, modifying the specified options with the specified values.
%
% options = cplexoptimset(oldopts,newopts) combines an existing options
% structure, oldopts, with a new options structure, newopts. Any options in
% newopts with nonempty values overwrite the corresponding old options in
% oldopts.
%
% Options corresponding to the MATLAB Optimization Toolbox
% Diagnostics          'on' | {'off'}
% Display              'off' | 'iter' | 'final' | 'notify'
% MaxIter               refer to cplex.Param.simplex.limits.iterations
% Simplex               refer to cplex.Param.lpmethod
%                       refer to cplex.Param.qpmethod
%                       refer to cplex.Param.mip.strategy.startalgorithm
% BranchStrategy        refer to cplex.Param.mip.strategy.variableselect
% MaxNodes              refer to cplex.Param.mip.limits.nodes
% MaxTime               refer to cplex.Param.timelimit
% NodeDisplayInterval   refer to cplex.Param.mip.interval
% NodeSearchStrategy    refer to cplex.Param.mip.strategy.nodeselect
% TolFun                refer to cplex.Param.simplex.tolerances.optimality
% TolXInteger           refer to cplex.Param.mip.tolerances.integrality
% TolRLPFun             refer to cplex.Param.simplex.tolerances.optimality
%
% For other CPLEX options, refer to the CPLEX Parameters Reference
% Manual.
%
% Examples
% To turn on the optimizer output and set a node limit of 400:
% options = cplexoptimset('Diagnostics', 'on', 'MaxNodes', 400);
%
% To set a node limit of 400 and instruct CPLEX to use traditional branch
% and cut style search:
% opt = cplexoptimset('cplex'); 
% opt.mip.limits.nodes=400; 
% opt.mip.strategy.search=1;
% 
%  See also cplexoptimget
%

% ---------------------------------------------------------------------------
% File: cplexoptimset.m
% ---------------------------------------------------------------------------
% Licensed Materials - Property of IBM
% 5725-A06 5725-A29 5724-Y48 5724-Y49 5724-Y54 5724-Y55
% Copyright IBM Corporation 2008, 2010. All Rights Reserved.
%
% US Government Users Restricted Rights - Use, duplication or
% disclosure restricted by GSA ADP Schedule Contract with
% IBM Corp.
% ---------------------------------------------------------------------------
