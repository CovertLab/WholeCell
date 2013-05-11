function [x,resnorm,residual,exitflag,output,lambda]=cplexlsqnonneglin(C,d,Aineq,bineq,Aeq,beq,x0,options)
%%
% Purpose
% Solve nonnegative least squares problems.
%
% Syntax
%    x = cplexlsqnonneglin(C,d)
%    x = cplexlsqnonneglin(C,d,Aineq,bineq)
%    x = cplexlsqnonneglin(C,d,Aineq,bineq,Aeq,beq)
%    x = cplexlsqnonneglin(C,d,Aineq,bineq,Aeq,beq,x0)
%    x = cplexlsqnonneglin(C,d,Aineq,bineq,Aeq,beq,x0,options)
%    x = cplexlsqnonneglin(problem)
%    [x,resnorm] = cplexlsqnonneglin(...)
%    [x,resnorm,residual] = cplexlsqnonneglin(...)
%    [x,resnorm,residual,exitflag] = cplexlsqnonneglin(...)
%    [x,resnorm,residual,exitflag,output] = cplexlsqnonneglin(...)
%    [x,resnorm,residual,exitflag,output,lambda] = cplexlsqnonneglin(...)
%
% Description
% Finds the minimum of a problem specified by
%    min      norm(C*x-d)^2
%    st.      Aineq*x <= bineq
%             Aeq*x    = beq
%             x >= 0
%
% d, bineq, beq are column vectors.
% C, Aineq, Aeq are matrices.
%
% x = cplexlsqnonneglin(C,d) solves the least squares problem min
% norm(C*x-d)^2 such that x >= 0.
%
% x = cplexlsqnonneglin(C,d,Aineq,bineq) solves the preceding problem 
% with the additional inequality constraints Aineq*x <= bineq.
%
% x = cplexlsqnonneglin(C,d,Aineq,bineq,Aeq,beq) solves the preceding
% problem with the additional equality constraints Aeq*x = beq. If no 
% inequalities exist, set Aineq=[] and bineq=[]. 
%
% x = cplexlsqnonneglin(C,d,Aineq,bineq,Aeq,beq,x0) sets the starting point 
% for the algorithm to x0.  If no equalities exist, set Aeq=[] and beq=[].
%
% x = cplexlsqnonneglin(C,d,Aineq,bineq,Aeq,beq,x0,options) minimizes with 
% the default optimization options replaced by values in the structure 
% options, which can be created using the function cplexoptimset. 
% If you do not wish to give an initial point, set x0=[].
%
% x = cplexlsqnoneglin(problem) where problem is a structure.
%
% [x,resnorm] = cplexlsqnonneglin(...) returns the value of the objective
% function at the solution x: resnorm = norm(C*x-d)^2.
%
% [x,resnorm,residual] = cplexlsqnonneglin(...) returns the residual at the
% solution: C*x-d.
%
% [x,resnorm,residual,exitflag] = cplexlsqnonneglin(...) returns a value
% exitflag that describes the exit condition of cplexlsqnonneglin.
%
% [x,resnorm,residual,exitflag,output] = cplexlsqnonneglin(...) returns a
% structure output that contains information about the optimization.
%
% [x,resnorm,residual,exitflag,output,lambda] = cplexlsqnonneglin(...)
% returns a structure lambda whose fields contain the Lagrange multipliers
% at the solution x.
%
% Input Arguments
% problem   Structure containing the following fields:
%           C 	    Double matrix for objective function
%           d       Double column vector for objective function
%           Aineq   Double matrix for linear inequality constraints
%           bineq   Double column vector for linear inequality
%                   constraints
%           Aeq     Double matrix for linear equality constraints
%           beq     Double column vector for linear equality constraints
%           x0      Double column vector for initial point of x
%           options Options structure created with cplexoptimset
%
% Output Arguments
% x         Solution found by the optimization function. If exitflag > 0,
%           then x is a solution; otherwise, x is the value of the
%           optimization routine when it terminated prematurely.
% resnorm   Value of the objective function at the solution x
% residual  Residual at the solution
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
% File: cplexlsqnonneglin.m
% ---------------------------------------------------------------------------
% Licensed Materials - Property of IBM
% 5725-A06 5725-A29 5724-Y48 5724-Y49 5724-Y54 5724-Y55
% Copyright IBM Corporation 2008, 2010. All Rights Reserved.
%
% US Government Users Restricted Rights - Use, duplication or
% disclosure restricted by GSA ADP Schedule Contract with
% IBM Corp.
% ---------------------------------------------------------------------------
