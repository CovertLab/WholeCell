function [x,resnorm,residual,exitflag,output]=cplexlsqmiqcp(C,d,Aineq,bineq,Aeq,beq,l,Q,r,sostype,sosind,soswt,lb,ub,ctype,x0,options)
%%
% Purpose
% Solve quadratically constrained mixed integer least squares problems.
%
% Syntax
%    x = cplexlsqmiqcp(C,d,Aineq,bineq)
%    x = cplexlsqmiqcp(C,d,Aineq,bineq,Aeq,beq)
%    x = cplexlsqmiqcp(C,d,Aineq,bineq,Aeq,beq,l,Q,r)
%    x = cplexlsqmiqcp(C,d,Aineq,bineq,Aeq,beq,l,Q,r,sostype,sosind,soswt)
%    x = cplexlsqmiqcp(C,d,Aineq,bineq,Aeq,beq,l,Q,r,sostype,sosind,soswt,
%        lb,ub)
%    x = cplexlsqmiqcp(C,d,Aineq,bineq,Aeq,beq,l,Q,r,sostype,sosind,soswt,
%        lb,ub,ctype)
%    x = cplexlsqmiqcp(C,d,Aineq,bineq,Aeq,beq,l,Q,r,sostype,sosind,soswt,
%        lb,ub,ctype,x0)
%    x = cplexlsqmiqcp(C,d,Aineq,bineq,Aeq,beq,l,Q,r,sostype,sosind,soswt,
%        lb,ub,ctype,x0,options)
%    x = cplexlsqmiqcp(problem)
%    [x,resnorm] = cplexlsqmiqcp(...)
%    [x,resnorm,residual] = cplexlsqmiqcp(...)
%    [x,resnorm,residual,exitflag] = cplexlsqmiqcp(...)
%    [x,resnorm,residual,exitflag,output] = cplexlsqmiqcp(...)
%
% Description
% Finds the minimum of a problem specified by
%    min      norm(C*x-d)^2
%    st.      Aineq*x      <= bineq
%             Aeq*x         = beq
%             l*x + x'*Q*x <= r
%             lb <= x <= ub
%             x belongs to BICSN
%
% d, bineq, beq, l, r, lb, and ub are column vectors.
% C, Aineq, Aeq and Q are matrices.
% x is a BICSN vector -- that is, its individual entries are each required
% to be binary, general integer, continuous, semi-continuous or
% semi-integer.
%
% x = cplexlsqmiqcp(C,d,Aineq,bineq) solves the constrained mixed integer
% least squares problem min norm(C*x-d)^2 such that Aineq*x <= bineq.
%
% x = cplexlsqmiqcp(C,d,Aineq,bineq,Aeq,beq) solves the preceding problem
% with the additional equality constraints Aeq*x = beq. If no inequalities
% exist, set Aineq=[] and bineq=[].
%
% x = cplexlsqmiqcp(C,d,Aineq,bineq,Aeq,beq,l,Q,r) solves the preceding
% problem while additionally satisfying the quadratic inequality
% constraints l*x + x'*Q*x <= r. If no equalities exist, set Aeq=[] 
% and beq=[].
%
% x = cplexlsqmiqcp(C,d,Aineq,bineq,Aeq,beq,l,Q,r,sostype,sosind,soswt)
% solves the preceding problem with the additional requirement that the SOS
% constraints are satisfied. If no quadratic inequalities exist, set l=[],
% Q=[] and r=[].
%
% x =
% cplexlsqmiqcp(C,d,Aineq,bineq,Aeq,beq,l,Q,r,sostype,sosind,soswt,lb,ub)
% defines a set of lower and upper bounds on the design variables, x, so
% that the solution is always in the range lb <= x <= ub. If no SOS
% constraints exist, set sostype=[],sosind=[] and soswt=[].
%
% x =
% cplexlsqmiqcp(C,d,Aineq,bineq,Aeq,beq,l,Q,r,sostype,sosind,soswt,lb,ub,ct
% ype) defines the types for each of the design variables. If no bounds
% exist, set lb=[] and ub=[].
%
% x =
% cplexlsqmiqcp(C,d,Aineq,bineq,Aeq,beq,l,Q,r,sostype,sosind,soswt,lb,ub,ct
% ype,x0) sets the starting point for the algorithm to x0. If all design
% variables are continuous, set ctype=[].
%
% x =
% cplexlsqmiqcp(C,d,Aineq,bineq,Aeq,beq,l,Q,r,sostype,sosind,soswt,lb,ub,ct
% ype,x0,options) minimizes with the default optimization options replaced
% by values in the structure options, which can be created using the
% function cplexoptimset. If you do not wish to give an initial point, set
% x0=[].
%
% x = cplexlsqmiqcp(problem) where problem is a structure.
%
% [x,resnorm] = cplexlsqmiqcp(...) returns the value of the objective
% function at the solution x: resnorm = norm(C*x-d)^2.
%
% [x,resnorm,residual] = cplexlsqmiqcp(...) returns the residual at the
% solution: C*x-d.
%
% [x,resnorm,residual,exitflag] = cplexlsqmiqcp(...) returns a value
% exitflag that describes the exit condition of cplexlsqmiqcp.
%
% [x,resnorm,residual,exitflag,output] = cplexlsqmiqcp(...) returns a
% structure output that contains information about the optimization.
%
% Input Arguments
% C         Double matrix for objective function
% d         Double column vector for objective function
% Aineq     Double matrix for linear inequality constraints
% bineq     Double column vector for linear inequality constraints
% Aeq       Double matrix for linear equality constraints
% beq       Double column vector for linear equality constraints
% l         Double column vector or matrix
%           Linear part of quadratic constraints
% Q         Double matrix or double matrix cell for quadratic
%           constraints
% r         Double or double row vector
%           Rhs of quadratic inequality constraints
% sostype   String with possible char values  '1', '2'
% sosind    Double column vector or column vector cell of indices for
%           the SOSs to be added
% soswt     Double column vector or column vector cell of weights for
%                     the SOSs to be added
% lb        Double column vector of lower bounds
% ub        Double column vector of upper bounds
% ctype     String with possible char values 'B','I','C','S','N'
%           ctype(j) to 'B', 'I','C', 'S', or 'N' to indicate
%           that x(j) should be binary, general integer,
%           continuous, semi-continuous or semi-integer (respectively).
% x0        Double column vector for initial point of x
% options   Options structure created with cplexoptimset
%
% problem   Structure containing the following fields:
%           C         Double matrix for objective function
%           d         Double column vector for objective function
%           Aineq     Double matrix for linear inequality constraints
%           bineq     Double column vector for linear inequality constraints
%           Aeq       Double matrix for linear equality constraints
%           beq       Double column vector for linear equality constraints
%      	    qc 	      Struct vector
%  	        qc(i).a 	 Double column vector for linear part of the 
%                        quadratic constraint
%  	        qc(i).rhs    Double for righthand side for quadratic constraint
%  	        qc(i).Q 	 Double matrix for the quadratic part of the
%                        quadratic constraint
%        	sos 	  Struct vector representing the SOSs
%  	        sos(i).type  String with possible char values  '1', '2'
%       	sos(i).ind   Double column vector of indices for the SOSs to
%                        be added
%  	        sos(i).wt    Double column vector of weights for the SOSs to 	
%                        be added
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
% File: cplexlsqmiqcp.m
% ---------------------------------------------------------------------------
% Licensed Materials - Property of IBM
% 5725-A06 5725-A29 5724-Y48 5724-Y49 5724-Y54 5724-Y55
% Copyright IBM Corporation 2008, 2010. All Rights Reserved.
%
% US Government Users Restricted Rights - Use, duplication or
% disclosure restricted by GSA ADP Schedule Contract with
% IBM Corp.
% ---------------------------------------------------------------------------
