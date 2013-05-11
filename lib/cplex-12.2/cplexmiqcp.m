function [x,fval,exitflag,output]=cplexmiqcp(H, f, Aineq, bineq, Aeq, beq, l, Q, r, sostype, sosind, soswt, lb, ub, ctype, x0, options)
%%
% Purpose
% Solve quadratically constrained linear or quadratic integer programming
% problems.
%
% Syntax
%    x = cplexmiqcp(H,f,Aineq,bineq)
%    x = cplexmiqcp(H,f,Aineq,bineq,Aeq,beq)
%    x = cplexmiqcp(H,f,Aineq,bineq,Aeq,beq,l,Q,r)
%    x = cplexmiqcp(H,f,Aineq,bineq,Aeq,beq,l,Q,r,sostype,sosind,soswt)
%    x = cplexmiqcp(H,f,Aineq,bineq,Aeq,beq,l,Q,r,sostype,sosind,soswt,lb,
%        ub)
%    x = cplexmiqcp(H,f,Aineq,bineq,Aeq,beq,l,Q,r,sostype,sosind,soswt,lb,
%        ub,ctype,x0)
%    x = cplexmiqcp(H,f,Aineq,bineq,Aeq,beq,l,Q,r,sostype,sosind,soswt,lb,
%        ub,ctype,x0,options)
%    x = cplexmiqcp(problem)
%    [x,fval] = cplexmiqcp(...)
%    [x,fval,exitflag] = cplexmiqcp(...)
%    [x,fval,exitflag,output] = cplexmiqcp(...)
%    [x,fval,exitflag,output,lambda] = cplexmiqcp(...)
%
% Description
% Finds the minimum of a problem specified by
%    min      0.5*x'*H*x+f*x or f*x
%    st.      Aineq*x      <= bineq
%             Aeq*x         = beq
%             l*x + x'*Q*x <= r
%             lb <= x <= ub
%             x belongs to BICSN
%
% f, bineq, beq, l, r, lb and ub are column vectors.
% H, Aineq, Aeq, and Q are matrices.
% x is a BICSN vector -- that is, its individual entries are each required
% to be binary, general integer, continuous, semi-continuous or 
% semi-integer.
%
% x = cplexmiqcp(H,f,Aineq,bineq) solves the mixed integer programming
% problem min 1/2*x'*H*x + f*x subject to Aineq*x <= bineq. If no
% quadratic objective term exists, set H=[].
%
% x = cplexmiqcp(H,f,Aineq,bineq,Aeq,beq) solves the preceding problem
% while additionally satisfying the equality constraints Aeq*x = beq. If no
% inequalities exist, set Aineq=[] and bineq=[].
%
% x = cplexmiqcp(H,f,Aineq,bineq,Aeq,beq,l,Q,r) solves the preceding
% problem while additionally satisfying the quadratic inequality
% constraints l*x + x'*Q*x <= r. If no equalities exist, set Aeq=[] and beq=[].
%
% x = cplexmiqcp(f,Aineq,bineq,Aeq,beq,l,Q,r,sostype,sosind,soswt) solves
% the preceding problem with the additional requirement that the SOS
% constraints are satisfied. If no quadratic inequalities exist, set l=[],
% Q=[] and r=[].
%
% x = cplexmiqcp(H,f,Aineq,bineq,Aeq,beq,l,Q,r,sostype,sosind,soswt,lb,ub)
% defines a set of lower and upper bounds on the design variables, x, so
% that the solution is in the range lb <= x <= ub. If no SOS constraints
% exist, set sostype=[],sosind=[] and soswt=[].
%
% x =
% cplexmiqcp(H,f,Aineq,bineq,Aeq,beq,l,Q,r,sostype,sosind,soswt,lb,ub,ctype
% ) defines the types for each of the design variables. If no bounds exist,
% set lb=[] and ub=[].
%
% x =
% cplexmiqcp(H,f,Aineq,bineq,Aeq,beq,l,Q,r,sostype,sosind,soswt,lb,ub,x0)
% sets the starting point to x0. If all design variables are continuous,
% set ctype=[].
%
% x = cplexmiqcp(H,f,Aineq,bineq,Aeq,beq,l,Q,r,lb,ub,ctype,x0,options) minimizes
% with the optimization options specified in the structure options, which
% can be created using the function cplexoptimset If you do not wish to
% give an initial point, set x0=[].
%
% x = cplexmiqcp(problem) where problem is a structure.
%
% [x,fval] = cplexmiqcp(...) returns the value of the objective function at
% the solution x: fval = 0.5*x'*H*x + f*x.
%
% [x,fval,exitflag] = cplexmiqcp(...) returns a value exitflag that
% describes the exit condition of cplexmiqcp.
%
% [x,fval,exitflag,output] = cplexmiqcp(...) returns a structure output
% that contains information about the optimization.
%
% Input Arguments
% H         Double matrix for objective function
% f         Double column vector for objective function
% Aineq     Double matrix for linear inequality constraints
% bineq     Double column vector for linear inequality constraints
% Aeq       Double matrix for linear equality constraints
% beq       Double column vector for linear equality constraints
% l         Double column vector or matrix
%           Linear part of quadratic constraints
% Q         Double matrix or double matrix cell for quadratic
%           constraints
% r         Double or double row vector
%           Righthand side of quadratic inequality constraints
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
%           (respectively).
% x0        Double column vector for initial point for x
% options   Options structure created with cplexoptimset
%
% problem   Structure containing the following fields:
%           H         Double matrix for objective function
%           f         Double column vector for objective function
%           Aineq     Double matrix for linear inequality constraints
%           bineq     Double column vector for linear inequality constraints
%           Aeq       Double matrix for linear equality constraints
%           beq       Double column vector for linear equality constraints
%      	    qc 	      Struct vector
%  	        qc(i).a 	 Double column vector for linear part of the 
%                        quadratic constraint
%  	        qc(i).rhs    Double for righthand side for quadratic constraint
%  	        qc(i).Q 	 Double matrix for quadratic part of the 
%                        quadratic constraint
%        	sos 	  Struct vector representing the SOSs
%  	        sos(i).type  String with possible char values  '1', '2'
%       	sos(i).ind   Double column vector of indices for the SOSs 	
%                        to be added
%  	        sos(i).wt    Double column vector of weights for the SOSs	
%                        to be added
%           lb        Double column vector of lower bounds
%           ub        Double column vector of upper bounds
%           ctype     String with possible char values 'B','I','C','S','N'
%                     ctype(j) to 'B', 'I','C', 'S', or 'N' to indicate
%                     that x(j) should be binary, general integer,
%                     continuous, semi-continuous or semi-integer
%                     (respectively).
%           x0        Double column vector for initial point for x
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
% File: cplexmiqcp.m
% ---------------------------------------------------------------------------
% Licensed Materials - Property of IBM
% 5725-A06 5725-A29 5724-Y48 5724-Y49 5724-Y54 5724-Y55
% Copyright IBM Corporation 2008, 2010. All Rights Reserved.
%
% US Government Users Restricted Rights - Use, duplication or
% disclosure restricted by GSA ADP Schedule Contract with
% IBM Corp.
% ---------------------------------------------------------------------------
