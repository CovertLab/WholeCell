function [x, fval, exitflag, output] = cplexbilp(f, Aineq, bineq, Aeq, beq, x0, options)
%%
% Purpose
% Solve binary integer programming problems.
%
% Syntax
%    x = cplexbilp(f)
%    x = cplexbilp(f, Aineq, bineq)
%    x = cplexbilp(f, Aineq, bineq, Aeq, beq)
%    x = cplexbilp(f, Aineq, bineq, Aeq, beq, x0)
%    x = cplexbilp(f, Aineq, bineq, Aeq, beq, x0, options)
%    x = cplexbilp(problem)
%    [x, fval] = cplexbilp(...)
%    [x, fval, exitflag] = cplexbilp(...)
%    [x, fval, exitflag, output] = cplexbilp(...)
%
% Description
% Finds the minimum of a problem specified by
%    min      f*x
%    st.      Aineq*x <= bineq
%             Aeq*x    = beq
%             x binary
%
% f, bineq, and beq are column vectors.
% Aineq and Aeq are matrices.
% x is required to be a binary integer vector--that is, its entries can
% only take on the values 0 or 1.
%
% x = cplexbilp(f) solves the binary integer programming problem min f*x.
%
% x = cplexbilp(f, Aineq, bineq) solves the binary integer programming
% problem min f*x such that Aineq*x <= bineq.
%
% x = cplexbilp(f, Aineq, bineq, Aeq, beq) solves the preceding problem
% with the additional equality constraints Aeq*x = beq. If no inequalities
% exist, set Aineq=[] and bineq=[].
%
% x = cplexbilp(f, Aineq, bineq, Aeq, beq, x0) sets the starting point for
% the algorithm to x0. If no equalities exist, set Aeq=[] and beq=[].
%
% x = cplexbilp(f, Aineq, bineq, Aeq, beq, x0, options) minimizes with the
% default optimization options replaced by values in the structure options,
% which can be created using the function cplexoptimset. If you do not wish
% to give an initial point, set x0=[].
%
% x = cplexbilp(problem) where problem is a structure.
%
% [x,fval] = cplexbilp(...) returns the value of the objective function at
% the solution x: fval = f*x.
%
% [x,fval,exitflag] = cplexbilp(...) returns a value exitflag that
% describes the exit condition of cplexbilp.
%
% [x,fval,exitflag,output] = cplexbilp(...) returns a structure output that
% contains information about the optimization.
%
% Input Arguments
% problem   Structure containing the following fields:
%           f         Double column vector for objective function
%           Aineq     Double matrix for linear inequality constraints
%           bineq     Double column vector for linear inequality 
%                     constraints
%           Aeq       Double matrix for linear equality constraints
%           beq       Double column vector for linear equality constraints
%           x0        Double column vector of initial point of x
%           options   Options structure created with cplexoptimset
%
% Output Arguments
% x         Solution found by the optimization function. If exitflag > 0,
%           then x is a solution; otherwise, x is the value of the
%           optimization routine when it terminated prematurely.
% fval      Value of the objective function at the solution x
% exitflag  Integer identifying the reason the optimization algorithm
%           terminated
% output    Structure containing information about the optimization. The
%           fields of the structure are:
%           iterations         Number of iterations
%           algorithm          Optimization algorithm used
%           message            Exit message
%           time               Execution time of the algorithm
%           cplexstatus        Status code of the solution
%           cplexstatusstring  Status string of the solution
%
%
%  See also cplexoptimset
%

% ---------------------------------------------------------------------------
% File: cplexbilp.m
% ---------------------------------------------------------------------------
% Licensed Materials - Property of IBM
% 5725-A06 5725-A29 5724-Y48 5724-Y49 5724-Y54 5724-Y55
% Copyright IBM Corporation 2008, 2010. All Rights Reserved.
%
% US Government Users Restricted Rights - Use, duplication or
% disclosure restricted by GSA ADP Schedule Contract with
% IBM Corp.
% ---------------------------------------------------------------------------
