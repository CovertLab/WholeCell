function [x,resnorm,residual,exitflag,output] = cplexlsqmilp (C,d,Aineq,bineq,Aeq,beq,sostype,sosind,soswt,lb,ub,ctype,x0,options)
%%
% Purpose
% Solve constrained mixed integer least squares problems.
%
% Syntax
%    x = cplexlsqmilp(C,d,Aineq,bineq)
%    x = cplexlsqmilp(C,d,Aineq,bineq,Aeq,beq)
%    x = cplexlsqmilp(C,d,Aineq,bineq,Aeq,beq,sostype,sosind,soswt)
%    x = cplexlsqmilp(C,d,Aineq,bineq,Aeq,beq,sostype,sosind,soswt,lb,ub)
%    x = cplexlsqmilp(C,d,Aineq,bineq,Aeq,beq,sostype,sosind,soswt,lb,ub,
%        ctype)
%    x = cplexlsqmilp(C,d,Aineq,bineq,Aeq,beq,sostype,sosind,soswt,lb,ub,
%        ctype,x0)
%    x = cplexlsqmilp(C,d,Aineq,bineq,Aeq,beq,sostype,sosind,soswt,lb,ub,
%        ctype,x0,options)
%    x = cplexlsqmilp(problem)
%    [x,resnorm] = cplexlsqmilp(...)
%    [x,resnorm,residual] = cplexlsqmilp(...)
%    [x,resnorm,residual,exitflag] = cplexlsqmilp(...)
%    [x,resnorm,residual,exitflag,output] = cplexlsqmilp(...)
%
% Description
% Finds the minimum of a problem specified by
%    min      norm(C*x-d)^2
%    st.      Aineq*x <= bineq
%             Aeq*x    = beq
%             lb <= x <= ub
%             x belongs to BICSN
%
% d, bineq, beq, lb, and ub are column vectors.
% C, Aineq, and Aeq are matrices.
% x is a BICSN vector -- that is, its individual entries are each required
% to be binary, general integer, continuous, semi-continuous or
% semi-integer.
%
% x = cplexlsqmilp(C,d,Aineq,bineq) solves the constrained mixed integer
% least squares problem min norm(C*x-d)^2 such that Aineq*x <= bineq.
%
% x = cplexlsqmilp(C,d,Aineq,bineq,Aeq,beq) solves the preceding problem
% with the additional equality constraints Aeq*x = beq. If no inequalities
% exist, set Aineq=[] and bineq=[].
%
% x = cplexlsqmilp(C,d,Aineq,bineq,Aeq,beq,sostype,sosind,soswt) solves the
% preceding problem with the additional requirement that the SOS
% constraints are satisfied. If no equalities exist, set Aeq=[] and beq=[].
%
% x = cplexlsqmilp(C,d,Aineq,bineq,Aeq,beq,sostype,sosind,soswt,lb,ub)
% defines a set of lower and upper bounds on the design variables, x, so
% that the solution is always in the range lb <= x <= ub. If no SOS
% constraints exist, set sostype=[],sosind=[] and soswt=[].
%
% x =
% cplexlsqmilp(C,d,Aineq,bineq,Aeq,beq,sostype,sosind,soswt,lb,ub,ctype)
% defines the types for each of the design variables. If no bounds exist,
% set lb=[] and ub=[].
%
% x =
% cplexlsqmilp(C,d,Aineq,bineq,Aeq,beq,sostype,sosind,soswt,lb,ub,ctype,x0)
% sets the starting point for the algorithm to x0. If all design variables
% are continuous, set ctype=[].
%
% x =
% cplexlsqmilp(C,d,Aineq,bineq,Aeq,beq,sostype,sosind,soswt,lb,ub,ctype,x0,
% options) minimizes with the default optimization options replaced by
% values in the structure options, which can be created using the function
% cplexoptimset. If you do not wish to give an initial point, set x0=[].
%
% x = cplexlsqmilp(problem) where problem is a structure.
%
% [x,resnorm] = cplexlsqmilp(...) returns the value of the objective
% function at the solution x: resnorm = norm(C*x-d)^2.
%
% [x,resnorm,residual] = cplexlsqmilp(...) returns the residual at the
% solution: C*x-d.
%
% [x,resnorm,residual,exitflag] = cplexlsqmilp(...) returns a value
% exitflag that describes the exit condition of cplexlsqmilp.
%
% [x,resnorm,residual,exitflag,output] = cplexlsqmilp(...) returns a
% structure output that contains information about the optimization.
%
% Input Arguments
% C         Double matrix for objective function
% d         Double column vector for objective function
% Aineq     Double matrix for linear inequality constraints
% bineq     Double column vector for linear inequality constraints
% Aeq       Double matrix for linear equality constraints
% beq       Double column vector for linear equality constraints
% sostype   String with possible char values  '1', '2'
% sosind    Double column vector or column vector cell of indices for
%           the SOSs to be added
% soswt     Double column vector or column vector cell of weights for
%           the SOSs to be added
% lb        Double column vector of lower bounds
% ub        Double column vector of upper bounds
% ctype     String with possible char values 'B','I','C','S','N'
%           ctype(j) to 'B', 'I','C', 'S', or 'N' to indicate
%           that x(j) should be binary, general integer,
%           continuous, semi-continuous or semi-integer
%          (respectively).
% x0        Double column vector for initial point of x
% options   Options structure created with cplexoptimset
%
% problem   Structure containing the following fields:
%           C         Double matrix for objective function
%           d         Double column vector for objective function
%           Aineq     Double matrix for linear inequality constraints
%           bineq     Double column vector for linear inequality 
%                     constraints
%           Aeq       Double matrix for linear equality constraints
%           beq       Double column vector for linear equality constraints
%        	sos 	  Struct vector representing the SOSs
%  	        sos(i).type    String with possible char values  '1', '2'
%       	sos(i).ind 	   Double column vector of indices for the SOSs to
%                          be added
%  	        sos(i).wt 	   Double column vector of weights for the SOSs to	
%                          be added
%           lb        Double column vector of lower bounds
%           ub        Double column vector of upper bounds
%           ctype     String with possible char values 'B','I','C','S','N'
%                     ctype(j) to 'B', 'I','C', 'S', or 'N' to indicate
%                     that x(j) should be binary, general integer,
%                     continuous, semi-continuous or semi-integer
%                     (respectively).
%           x0        Double column vector for initial point of x
%           options   Options structure created with cplexoptimset
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
%
%
%  See also cplexoptimset
%

% ---------------------------------------------------------------------------
% File: cplexlsqmilp.m
% ---------------------------------------------------------------------------
% Licensed Materials - Property of IBM
% 5725-A06 5725-A29 5724-Y48 5724-Y49 5724-Y54 5724-Y55
% Copyright IBM Corporation 2008, 2010. All Rights Reserved.
%
% US Government Users Restricted Rights - Use, duplication or
% disclosure restricted by GSA ADP Schedule Contract with
% IBM Corp.
% ---------------------------------------------------------------------------
