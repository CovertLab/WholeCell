function [x,fval,exitflag,output,lambda]=cplexqp(H,f,Aineq,bineq,Aeq,beq,lb,ub,x0,options)
%%
% Purpose
% Solve quadratic programming problems.
%
% Syntax
%   x = cplexqp(H,f,Aineq,bineq)
%   x = cplexqp(H,f,Aineq,bineq,Aeq,beq)
%   x = cplexqp(H,f,Aineq,bineq,Aeq,beq,lb,ub)
%   x = cplexqp(H,f,Aineq,bineq,Aeq,beq,lb,ub,x0)
%   x = cplexqp(H,f,Aineq,bineq,Aeq,beq,lb,ub,x0,options)
%   x = cplexqp(problem)
%   [x,fval] = cplexqp(...)
%   [x,fval,exitflag] = cplexqp(...)
%   [x,fval,exitflag,output] = cplexqp(...)
%   [x,fval,exitflag,output,lambda] = cplexqp(...)
%
% Description
% Finds the minimum of a problem specified by
%    min      0.5*x'*H*x+f*x or f*x
%    st.      Aineq*x <= bineq
%             Aeq*x    = beq
%             lb <= x <= ub
%
% f, bineq, beq, lb, and ub are column vectors.
% H, Aineq, and Aeq are matrices.
%
% x = cplexqp(H,f,Aineq,bineq) solves the quadratic programming problem min
% 1/2*x'*H*x + f*x subject to Aineq*x <= bineq.
%
% x = cplexqp(H,f,Aineq,bineq,Aeq,beq) solves the preceding problem with
% the additional equality constraints Aeq*x = beq. If no inequalities
% exist, set Aineq=[] and bineq=[].
%
% x = cplexqp(H,f,Aineq,bineq,Aeq,beq,lb,ub) defines a set of lower and
% upper bounds on the design variables, x, so that the solution is in the
% range lb <= x <= ub. If no equalities exist, set Aeq=[] and beq=[].
%
% x = cplexqp(H,f,Aineq,bineq,Aeq,beq,lb,ub,x0) sets the starting point to
% x0. If no bounds exist, set lb=[] and ub=[].
%
% x = cplexqp(H,f,Aineq,bineq,Aeq,beq,lb,ub,x0,options) minimizes with the
% default optimization options replaced by values in the structure options,
% which can be created using the function cplexoptimset. If you do not wish
% to give an initial point, set x0=[].
%
% x = cplexqp(problem) where problem is a structure.
%
% [x,fval] = cplexqp(...) returns the value of the objective function at
% the solution x: fval = 0.5*x'*H*x + f*x.
%
% [x,fval,exitflag] = cplexqp(...) returns a value exitflag that describes
% the exit condition of cplexqp.
%
% [x,fval,exitflag,output] = cplexqp(...) returns a structure output that
% contains information about the optimization.
%
% [x,fval,exitflag,output,lambda] = cplexqp(...) returns a structure
% lambda whose fields contain the Lagrange multipliers at the solution x.
%
% Input Arguments
% problem   Structure containing the following fields:
%           H         Double matrix for objective function
%           f         Double column vector for objective function
%           Aineq     Double matrix for linear inequality constraints
%           bineq     Double column vector for linear inequality 
%                     constraints
%           Aeq       Double matrix for linear equality constraints
%           beq       Double column vector for linear equality constraints
%           lb        Double column vector of lower bounds
%           ub        Double column vector of upper bounds
%           x0        Double column vector for initial point of x
%           options   Options structure created with cplexoptimset
%
% Output Arguments
% x         Solution found by the optimization function. If exitflag > 0,
%           then x is a solution; otherwise, x is the value of the
%           optimization routine when it terminated prematurely.
% fval      Value of the objective function at the solution x
% exitflag  Integer identifying the reason the optimization algorithm
%           terminated.
% output    Structure containing information about the optimization. The
%           fields of the structure are:
%           iterations         Number of iterations
%           algorithm          Optimization algorithm used
%           message            Exit message
%           time               Execution time of the algorithm
%           cplexstatus        Status code of the solution
%           cplexstatusstring  Status string of the solution
% lambda    Structure containing the Lagrange multipliers at the solution x
%           (separated by constraint type). The fields of the structure
%           are:
%           lower              Lower bounds lb
%           upper              Upper bounds ub
%           ineqlin            Linear inequalities
%           eqlin              Linear equalities
%
%
%  See also cplexoptimset
%

% ---------------------------------------------------------------------------
% File: cplexqp.m
% ---------------------------------------------------------------------------
% Licensed Materials - Property of IBM
% 5725-A06 5725-A29 5724-Y48 5724-Y49 5724-Y54 5724-Y55
% Copyright IBM Corporation 2008, 2010. All Rights Reserved.
%
% US Government Users Restricted Rights - Use, duplication or
% disclosure restricted by GSA ADP Schedule Contract with
% IBM Corp.
% ---------------------------------------------------------------------------
