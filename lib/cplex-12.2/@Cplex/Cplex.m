classdef Cplex < dynamicprops
   % Cplex The math programming solver.
   %   This class stores MP models and provides methods for its solution
   %   analysis and manipulation.
   %   Cplex properties:
   %     Model          The Model
   %     Param          The Parameters to solve the Model
   %     DisplayFunc    The display function handle of the Model
   %     Conflict       A dynamic property of the Cplex class that
   %                    represents the Conflict of the model,
   %                    it will be generated after refineConflict.
   %                    If the Model is MIP, refineMipStartConflict can
   %                    generate Conflict as well.
   %                    use cplex.findprop('Conflict').delete to remove
   %                    this property
   %     InfoCallback   A dynamic property, an informational callback is
   %                    a user-written routine that enables your
   %                    application to access
   %                    information about the current mixed integer
   %                    programming (MIP) optimization without sacrificing
   %                    performance and without interfering in the
   %                    search of the solution space.
   %                    use cplex.findprop('InfoCallback').delete to remove
   %                    this property
   %     Start          A dynamic property of the Cplex class that
   %                    represents the Start of the LP or QP model,
   %                    it will be generated after solving.
   %                    It also can be given before
   %                    solving to help CPLEX(R) get a Solution. It is
   %                    optional.
   %                    use cplex.findprop('Start').delete to remove
   %                    this property
   %     MipStart       A dynamic property of the Cplex class that
   %                    represents the MipStart of the MIP model;
   %                    it will be generated after solving
   %                    or populating. It also can be given before
   %                    solving to help CPLEX get a Solution. It is
   %                    optional.
   %                    use cplex.findprop('MipStart').delete to remove
   %                    this property
   %     Solution       A dynamic property of the Cplex class that
   %                    represents the solution of the model,
   %                    it will be generated after solving or feasopting.
   %                    If the Model is MIP, populate can generate
   %                    Solution as well.
   %                    use cplex.findprop('Solution').delete to remove
   %                    this property
   %   Cplex methods:
   %     Cplex          The constructor for Cplex objects.
   %     addCols        Adds columns to the problem object.
   %     addIndicators  Adds indicator constraints to the specified problem.
   %     addQCs         Adds a quadratic constraint or quadratic
   %                    constraint set to the problem object.
   %     addRows        Adds constraints to the problem object.
   %     addSOSs        Adds information about special ordered sets (SOS)
   %                    to a problem object of type MILP, MIQP, or MIQCP.
   %     delCols        Removes the specified columns from the problem
   %                    object.
   %     delRows        Removes the specified rows from the problem
   %                    object.
   %     feasOpt        Computes a minimum-cost relaxation of the righthand
   %                    side values of constraints or bounds on variables
   %                    in order to make an infeasible problem feasible.
   %     getChgParam    Returns a struct of parameters which are not set at
   %                    their default values.
   %     getProbType    Accesses the problem type.
   %     getVersion     Returns the version of CPLEX.
   %     populate       Generates multiple Solutions to a mixed integer
   %                    programming (MIP) problem.
   %     readBasis      Reads a basis from a BAS file and copies that basis
   %                    into the Start property of the Model.
   %     readMipStart   Reads a MST file and copies the information of all
   %                    the MIP starts contained in this file into the
   %                    problem object.
   %     readModel      Reads a CPLEX Model from a file and copies it into
   %                    the problem object.
   %     readParam      Reads parameter names and settings from the file
   %                    specified by filename and copies them into the
   %                    problem object.
   %     refineConflict Identifies a minimal conflict for the infeasibility
   %                    of the linear constraints and the variable bounds in
   %                    the current linear problem.
   %     refineMipStartConflict
   %                    Refines a conflict in order to determine why a
   %                    given MIP start is not feasible in linear
   %                    problem.
   %     setDefault     Resets all CPLEX parameters and settings to default
   %                    values.
   %     solve          Solves the current model in the problem object.
   %     terminate      Interrupts an ongoing solve() method invocation.
   %     tuneParam      Tunes the Parameters of the environment for improved
   %                    optimizer performance on the problem object.
   %     writeBasis     Writes the most current basis associated with the
   %                    problem object to a file.
   %     writeConflict  Writes a conflict file.
   %     writeMipStart  Writes a range of MIP starts of the problem
   %                    object to a file in MST format.
   %     writeModel     Writes the model of the invoking object to a
   %                    file.
   %     writeParam     Writes the parameter names and their current settings
   %                    into the file specified by name for all the CPLEX
   %                    parameters that are not currently set at their
   %                    default.
   %
   %   See also
   %   Properties:
   %      Model, Param, DisplayFunc
   %   Methods:
   %      Cplex
   %      getVersion     getProbType     getChgParam
   %      readBasis      writeBasis      readMipStart     writeMipStart
   %      readModel      writeModel      writeConflict
   %      readParam      writeParam      tuneParam        setDefault
   %      addCols        addRows         delRows          delCols
   %      addIndicators  addQCs          addSOSs
   %      solve          feasOpt         populate
   %      refineConflict terminate       refineMipStartConflict
   
% ---------------------------------------------------------------------------
% File: Cplex.m
% ---------------------------------------------------------------------------
% Licensed Materials - Property of IBM
% 5725-A06 5725-A29 5724-Y48 5724-Y49 5724-Y54 5724-Y55
% Copyright IBM Corporation 2008, 2010. All Rights Reserved.
%
% US Government Users Restricted Rights - Use, duplication or
% disclosure restricted by GSA ADP Schedule Contract with
% IBM Corp.
% ---------------------------------------------------------------------------
   properties
      % Model The Model used by the math programming solver.
      %   Model fields:
      %     name     String Name of the Model.
      %     sense    String that specifies whether the problem is a
      %              minimization or maximization problem.
      %     obj      Double column vector
      %              representing objective function coefficients.
      %     lb       Double column vector
      %              representing lower bound on each of the variables.
      %     ub       Double column vector
      %              representing upper bound on each of the variables.
      %     A        Constraint matrix.
      %     lhs      Double column vector
      %              representing lefthand side value for each constraint
      %              in the constraint matrix.
      %     rhs      Double column vector
      %              representing righthand side value for each constraint
      %              in the constraint matrix.
      %
      %   Optional fields of cplex.Model:
      %     objname   String representing the name of the objective.
      %     colname   Char matrix representing the names of the matrix
      %               columns or, equivalently, the variable names.
      %               Set colname by cplex.Model.colname(1,:) = 'mycolname'
      %     rowname   Char matrix representing the names of the matrix
      %               rows or, equivalently, the constraint names.
      %               Set colname by cplex.Model.rowname(1,:) = 'myrowname'
      %     ctype     String containing the type of each column in the
      %               constraint matrix, specifying whether a variable
      %               is continuous, integer, binary, semi-continuous or
      %               semi-integer. The possible char vales are
      %               'B','I','C','S','N'. Set ctype(j) to 'B', 'I','C',
      %               'S', or 'N' to indicate that x(j) should be binary,
      %               general integer, continuous, semi-continuous or
      %               semi-integer (respectively).
      %
      %     sos     Struct vector representing the SOSs.
      %     sos(i).name   (Optional) String representing the name of sos(i).
      %     sos(i).type   Type of sos(i), char '1' or '2'.
      %     sos(i).ind    Double column vector representing
      %                   the index of SOS to be added.
      %     sos(i).wt     Double column vector representing
      %                   the weights of the SOS to be added.
      %     Q      Quadratic objective matrix.
      %     qc     Struct vector representing the quadratic constraints.
      %     qc(i).name   (Optional) String representing the name of qc(i).
      %     qc(i).sense   Sense of quadratic constraint, char 'L' or 'G'.
      %     qc(i).a       Double column vector representing
      %                   the linear part of the quadratic constraint.
      %     qc(i).rhs     Righthand side term for the constraint.
      %     qc(i).Q       matrix representing the
      %                   Quadratic part of the quadratic constraint.
      %
      %     indicator Struct vector representing the indicator constraints.
      %     indicator(i).name       (Optional) String representing the name
      %                             of indicator(i).
      %     indicator(i).variable   Binary variable that acts as the
      %                             indicator for the constraint.
      %     indicator(i).complement Boolean value that specifies whether the
      %                             indicator variable is complemented.
      %     indicator(i).a      Double column representing
      %                         the linear portion of the indicator constraint.
      %     indicator(i).rhs    Righthand side value for the linear portion
      %                         of the indicator constraint.
      %     indicator(i).sense  char 'L','G' or 'E'.
      %     See also addIndicators, addSOSs, addQCs.
      Model;
      % Param Parameters for the math programming solver.
      %   The behavior of CPLEX(R) is controlled by a variety of parameters that
      %   can be accessed and modified by the user. Each parameter is a field
      %   of the Cplex.Param structure.
      %
      %   Example:
      %     Set a parameter:
      %       cplex.Param.lpmethod.Cur = 5
      %     Get a parameter:
      %       cplex.Param.lpmethod.Cur
      %
      %     Set a node limit of 400 and instruct CPLEX to use traditional 
      %     branch and cut style search:
      %       cplex.Param.mip.limits.nodes.Cur=400; 
      %       cplex.Param.mip.strategy.search.Cur=1;
      %
      %   For the usage of each parameter, please refer to the CPLEX
      %   Parameters Reference Manual
      %
      Param;
      % DisplayFunc   The Display Function of the Model
      %   The default value of this property is @disp, the function handle
      %   of the display function in MATLAB. With the default, all of the
      %   log information from CPLEX(R) will be displayed. If the
      %   Cplex.DisplayFunc property is set to empty, then the log
      %   information from CPLEX will not be displayed. In addition,
      %   users can write a custom DisplayFunc to control the output.
      DisplayFunc;
   end%end of properties
   properties (Hidden = true)
      Handle;
      Version;
   end%end of properties
   properties (Constant = true, Hidden = true)
      % Default is a struct contain the hierarchy of Params and Quality
      % And it is a Constant property in Cplex class, which mean it is a
      % property belong to Cplex class, not belong to any Cplex object,
      % then we can reuse Default in each Cplex object, which like the
      % static property in C++ class
      Default = cplexlink122 ([], 'getParamHierarchy');
   end%end of properties
      function out = display(cpx)
         % Overwrite  build-in display function
         % See also Cplex.
      end %end of function
      function out = terminate(cpx)
         % Cplex.terminate  Interrupts an ongoing solve() method invocation.
         %   Use this to implement stop buttons for a GUI.
         %
         %   See also Cplex.
      end %end of function
      function out = Cplex(modelname)
         % Cplex  The constructor for Cplex objects.  The constructor takes one
         %   optional argument to specify the problem name.
         %
         %     Example:
         %       cplex=Cplex('prob')
         %
         %   See also Cplex.
      end %end of function
      function out = getVersion(cpx)
         % Cplex.getVersion  Returns the version of CPLEX(R)
         %
         %   Example:
         %     cplex.getVersion()
         %
         %   See also Cplex.
      end %end of function
      function out = getChgParam(cpx)
         % Cplex.getChgParam  Returns a struct of parameters which are not 
         % set at their default values.
         %
         %   Example:
         %     cplex.getChgParam()
         %
         %   See also Cplex.
      end %end of function
      function out = readModel(cpx, filename)
         % Cplex.readModel  Reads a CPLEX(R) model from a file and copies it into
         %   the invoking object.
         %
         %   Example:
         %     cplex.readModel('myprob.lp')
         %
         %   For an example, refer to lpex2.m in the examples directory.
         %
         %   Parameters:
         %     filename    string
         %            The filename must end in one of these suffixes: .lp,
         %            .mps, .sav and .gz.
         %
         %   See also Cplex.
      end %end of function
      function out = writeModel(cpx, filename)
         % Cplex.writeModel  Writes the CPLEX(R) model of the invoking
         % object to a file
         %
         %   Example:
         %     cplex.writeModel('myprob.lp')
         %
         %   For an example, refer to lpex2.m in the examples directory.
         %
         %   Parameters:
         %     filename    string
         %            A file with one of the formats
         %              MPS     MPS format
         %              LP      CPLEX LP format with names modified to
         %                   conform to LP format
         %              SAV     Binary matrix and basis file 
         %              REW     MPS format with all names changed to
         %                   generic names
         %              RMP     MPS format, with all names changed to
         %                   generic names
         %              RLP     LP format, with all names changed to
         %                   generic names
         %
         %   See also Cplex.
      end %end of function
      function out = writeConflict(cpx, filename)
         % Cplex.writeConflict  Writes a Conflict file named filename.
         %
         %   Example:
         %     cplex.writeConflict('Conflict.lp')
         %
         %   Parameters:
         %     filename    string
         %            A file in the .lp format.
         %
         %   See also Cplex.
      end %end of function
      function out = readParam(cpx, filename)
         % Cplex.readParam   Reads parameters names and settings from the file
         %   specified by filename and copies them into the Cplex object.
         %
         %   This routine reads and copies files in the PRM format, as 
         %   created by writeParam. 
         %   The PRM format is documented in the CPLEX File Formats
         %   manual.
         %
         %   Example:
         %     cplex.readParam('myprob.prm')
         %
         %   Parameters:
         %     filename    string
         %            A file in the .prm format.
         %
         %   See also Cplex.
      end %end of function
      function out = writeParam(cpx, filename)
         % Cplex.writeParam   Writes the Parameter names and their current
         %   settingsinto the file specified by name for all the CPLEX(R)
         %   parameters that are not currently set at their default.
         %
         %   This routine reads and copies files in the PRM format, 
         %   as created by writeParam. The PRM format is documented in the 
         %   CPLEX File Formats manual.
         %
         %   Example:
         %     cplex.writeParam('myprob.prm')
         %
         %   Parameters:
         %     filename    string
         %            A file in the .prm format.
         %
         %   See also Cplex.
      end %end of function
      function out = readBasis(cpx, filename)
         % Cplex.readBasis   Reads a basis from a BAS file and copies that basis
         %   into a Cplex problem object.
         %
         %   The Parameter advance must be set to 1 (one), its default value,
         %   or 2(two)in order for the basis to be used for starting a
         %   subsequent optimization.
         %
         %   Example:
         %     cplex.readBasis('myprob.bas')
         %
         %   Parameters:
         %     filename    string
         %            A file in the .bas format.
         %
         %   See also Cplex.
      end %end of function
      function out = writeBasis(cpx, filename)
         % Cplex.writeBasis
         %   Writes the most current basis associated with a Cplex
         %   problem object to a file.
         %
         %   The file is saved in BAS format which corresponds to the industry
         %   standard MPS insert format for bases.
         %
         %   When writeBasis is invoked, the current basis is written 
         %   to a file.
         %   This routine does not remove the basis from the problem
         %   object.
         %
         %   Example:
         %     cplex.writeBasis('myprob.bas')
         %
         %   Parameters:
         %     filename    string
         %            A file in the .bas format.
         %
         %   See also Cplex.
      end %end of function
      function out = readMipStart(cpx, filename)
         % Cplex.readMipStart  Reads a MST file and copies the information of
         % all the MIP starts contained in this file into a Cplex problem
         % object.
         %
         %   The parameter advance must be set to 1 (one), its default value,
         %   or 2(two) in order for the MIP starts to be used.
         %
         %   Example:
         %     cplex.readMipStart('myprob.mst')
         %
         %   Parameters:
         %     filename    string
         %            A file in the .mst format.
         %
         %   See also Cplex.
      end %end of function
      function out = writeMipStart(cpx, filename)
         % Cplex.writeMipStart  Writes a range of MIP starts of a Cplex
         %   problem object to a file in MST format.
         %
         %   The MST format is an XML format and is documented in the
         %   stylesheet solution.xsl and schema solution.xsd in the
         %   include directory of the CPLEX(R) distribution. CPLEX
         %   File Formats also documents this format briefly.
         %
         %   Example:
         %     cplex.writeMipStart('myprob.mst')
         %
         %   Parameters:
         %     filename    string
         %            A file in the .mst format.
         %
         %   See also Cplex.
      end %end of function
      function out = solve(cpx)
         % Cplex.solve  Solves the current model in the invoking object.
         %
         % After a call to solve, the field Solution of the invoking
         % object will be populated as well as the field Start for an LP
         % or MipStart for a MIP.
         %
         %   Example:
         %     cplex.solve()
         %
         %   For an example, refer to lpex1.m in the examples directory.
         %
         %   See also Cplex.
      end %end of function
      function out = setDefault(cpx)
         % Cplex.setDefault  Resets all CPLEX(R) parameters and settings to
         %   default values
         %
         %   Example:
         %     cplex.setDefault()
         %
         %   See also Cplex.
      end %end of function
      function out = populate(cpx)
         % Cplex.populate  Generates multiple solutions to a mixed integer
         %   programming (MIP) problem.
         %
         % After a call to solve, the field Solution of the invoking
         % object will be populated as well as the field Start for an LP
         % or MipStart for a MIP.
         %
         %   Example:
         %     cplex.populate()
         %
         %   See also Cplex.
      end %end of function
      function out = feasOpt(cpx, preflhs, prefrhs, preflb, prefub)
         % Cplex.feasOpt  Computes a minimum-cost relaxation of the righthand
         %   side values of constraints or bounds on variables in order
         %   to make an infeasible problem feasible.
         %
         %   Example:
         %     cplex.feasOpt(preflhs, prefrhs, preflb, prefub)
         %     If one parameter is empty, set it to [].
         %     prefrhs=[]
         %     cplex.feasOpt(preflhs, prefrhs, preflb, prefub)
         %
         %   Parameters:
         %     preflhs    double column vector
         %                The length must be at least equal to the
         %                number of rows in the problem. An empty
         %                vector may be specified if no range values
         %                are allowed to be relaxed or none are
         %                present in the active problem.
         %                When not empty, the vector specifies the
         %                preference values that determine the cost
         %                of relaxing each range.
         %     prefrhs    double column vector
         %                The length must be at least equal to the
         %                number of rows in the problem. An empty
         %                vector may be specified if no rhs values
         %                are allowed to be relaxed.
         %                When not empty, the vector specifies the
         %                preference values that determine the cost
         %                of relaxing each constraint.
         %     preflb     double column vector
         %                The length must be at least equal to the
         %                number of columns in the problem. An empty
         %                vector may be passed if no lower bound of
         %                any variable is allowed to be relaxed.
         %                When not NULL, the vector specifies the
         %                preference values that determine the cost
         %                of relaxing each lower bound.
         %     prefub     double column vector
         %                The length must be at least equal to the
         %                number of columns in the problem. An empty
         %                vector may be passed if no upper bound of
         %                any variable is allowed to be relaxed.
         %                When not empty, the vector specifies the
         %                preference values that determine the cost
         %                of relaxing each upper bound.
         %
         %   Parameter FeasOptMode values:
         %   0   Minimize the sum of all required relaxations in first
         %       phase only; default
         %   1   Minimize the sum of all required relaxations in first
         %       phase and execute second phase to find optimum among
         %       minimal relaxations
         %   2   Minimize the number of constraints and bounds requiring
         %       relaxation in first phase only
         %   3   Minimize the number of constraints and bounds requiring
         %       relaxation in first phase and execute second phase to
         %       find optimum among minimal relaxations
         %   4   Minimize the sum of squares of required relaxations in
         %       first phase only
         %   5   Minimize the sum of squares
         %       of required relaxations in first phase and execute
         %       second phase to find optimum among minimal relaxations
         %
         %    It can minimize the weighted sum of the penalties for relaxations
         %    (denoted by SUM).
         %    It can minimize the weighted number of relaxed bounds and
         %    constraints (denoted by INF).
         %    It can minimize the weighted sum of the squared penalties of the
         %    relaxations (denoted by QUAD).
         %
         %   See also Cplex.
      end %end of function
      function out = refineMipStartConflict(cpx,index)
         % Cplex.refineMipStartConflict  Refines a Conflict in order to 
         %   determine why a given MIP start is not feasible in linear
         %   problem.
         %   In other words, this routine identifies a minimal conflict for the
         %   infeasibility of the linear constraints and bounds in a MIP start.
         %
         %   Example:
         %     cplex.refineMipStartConflict(index)
         %
         %   Parameters:
         %     index     Index of the MIP start among all the MIP starts
         %               associated with the problem.
         %   Returns:
         %     Cplex.Conflict  A struct vector with the fields:
         %           colind    Vector to receive the list of the
         %                     indices of the variables that
         %                     participate in the Conflict.
         %                     The length of the vector must not be
         %                     less than the number of columns in the
         %                     conflict.
         %                     If that number is not known, use the
         %                     number of columns in the problem
         %                     object instead.
         %         colbdstat   Vector to receive the Conflict
         %                     status of the columns.
         %                     Entry colbdstat(i) gives the status of
         %                     column colind(i).
         %                     The length of the vector must not be
         %                     less than the number of columns
         %                     in the Conflict.
         %                     If that number is not known, use the
         %                     number of columns in the problem
         %                     object instead.
         %          rowind     Vector to receive the list of the
         %                     indices of the constraints that
         %                     participate in the Conflict.
         %                     The length of the vector must not be
         %                     less than the number of rows in the
         %                     conflict.
         %                     If that number is not known, use the
         %                     total number of rows in the probem
         %                     object instead.
         %        rowbdstat    Vector to receive the Conflict
         %                     status of the rows.
         %                     Entry rowbdstat(i) gives the status of
         %                     row rowind(i).
         %                     The length of the vector must not be
         %                     less than the number of rows in the
         %                     conflict.
         %                     If that number is not known, use the
         %                     number of rows in the problem object
         %                     instead.
         %           status    Status of the Conflict
         %
         %   Conflict status values:
         %   Status     Meaning
         %    30    The problem appears to be feasible; no conflict
         %          is available.
         %    31    The conflict refiner found a minimal conflict.
         %    32    The conflict refiner concluded contradictory
         %          feasibility for the same set of constraints
         %          due to numeric problems. A conflict is
         %          available, but it is not minimal.
         %    33    The conflict refiner terminated because of a
         %          time limit. A conflict is available,
         %          but it is not minimal.
         %    34    The conflict refiner terminated because of an
         %          iteration limit. A conflict is available, but
         %          it is not minimal.
         %    35    The conflict refiner terminated because of a
         %          node limit. A conflict is available,
         %          but it is not minimal.
         %    36    The conflict refiner terminated because of an
         %          objective limit. A conflict is available,
         %          but it is not minimal.
         %    37    The conflict refiner terminated because of a
         %          memory limit. A conflict is available,
         %          but it is not minimal.
         %    38    The conflict refiner terminated because a user
         %          terminated the application. A conflict is
         %          available, but it is not minimal.
         %
         %   See also Cplex.
      end %end of function
      function out = refineConflict(cpx)
         % Cplex.refineConflict  Identifies a minimal Conflict for the
         %   infeasibility of the linear constraints and the variable bounds
         %   in the current
         %   linear problem.
         %
         %   Example:
         %     cplex.refineConflict()
         %
         %   Returns:
         %     Cplex.Conflict  A struct vector with the fields:
         %           colind    Vector to receive the list of the
         %                     indices of the variables that
         %                     participate in the Conflict.
         %                     The length of the vector must not be
         %                     less than the number of columns in the
         %                     Conflict.
         %                     If that number is not known, use the
         %                     number of columns in the problem
         %                     object instead.
         %         colbdstat   Vector to receive the Conflict
         %                     status of the columns.
         %                     Entry colbdstat[i] gives the status of
         %                     column colind[i].
         %                     The length of the vector must not be
         %                     less than the number of columns
         %                     in the Conflict.
         %                     If that number is not known, use the
         %                     number of columns in the problem
         %                     object instead.
         %          rowind     Vector to receive the list of the
         %                     indices of the constraints that
         %                     participate in the Conflict.
         %                     The length of the vector must not be
         %                     less than the number of rows in the
         %                     Conflict.
         %                     If that number is not known, use the
         %                     total number of rows in the problem
         %                     object instead.
         %         rowbdstat   Vector to receive the Conflict
         %                     status of the rows.
         %                     Entry rowbdstat[i] gives the status of
         %                     row rowind[i].
         %                     The length of the vector must not be
         %                     less than the number of rows in the
         %                     Conflict.
         %                     If that number is not known, use the
         %                     number of rows in the problem object
         %                     instead.
         %           status    Status of the Conflict.
         %
         %   Conflict status values:
         %   Status     Meaning
         %    30    The problem appears to be feasible; no conflict
         %          is available.
         %    31    The conflict refiner found a minimal conflict.
         %    32    The conflict refiner concluded contradictory
         %          feasibility for the same set of constraints
         %          due to numeric problems. A conflict is
         %          available,but it is not minimal.
         %    33    The conflict refiner terminated because of a
         %          time limit. A conflict is available,
         %          but it is not minimal.
         %    34    The conflict refiner terminated because of an
         %          iteration limit. A conflict is available, but
         %          it is not minimal.
         %    35    The conflict refiner terminated because of a
         %          node limit. A conflict is available,
         %          but it is not minimal.
         %    36    The conflict refiner terminated because of an
         %          objective limit. A conflict is available,
         %          but it is not minimal.
         %    37    The conflict refiner terminated because of a
         %          memory limit. A conflict is available,
         %          but it is not minimal.
         %    38    The conflict refiner terminated because a user
         %          terminated the application. A conflict is
         %          available, but it is not minimal.
         %
         %   See also Cplex.
      end %end of function
      function out = getProbType(cpx)
         % Cplex.getProbType  Accesses the problem type
         %
         %    Probtype     Meaning
         %     -1       Error: no problem or environment.
         %      0       Linear program; no quadratic data or ctype
         %              information stored.
         %      1       Problem with ctype information.
         %      3       Problem with ctype information, integer variables
         %              fixed.
         %      5       Problem with quadratic data stored.
         %      7       Problem with quadratic data and ctype information.
         %      8       Problem with quadratic data and ctype
         %              information, integer variables fixed.
         %     10       Problem with quadratic constraints.
         %     11       Problem with quadratic constraints and ctype
         %              information.
         %
         %   See also Cplex.
      end %end of function
      function out  = tuneParam(cpx,varargin)
         % Cplex.tuneParam  Tunes the Parameters of the environment for 
         % improved optimizer performance on the specified problem object.
         %
         %   Example:
         %    cplex.tuneParam();
         %    cplex.tuneParam({cplex.Param.mip.strategy.heuristicfreq,...
         %                    cplex.Param.mip.strategy.branch)
         %    cplex.tuneParam(cplex.Param.mip.cuts)
         %
         %   See also Cplex
      end %end of function
      function out = addRows(cpx, lhs, A, rhs, rowname)
         % Cplex.addRows  Adds constraints to a specified Cplex problem object.
         %
         %   Example:
         %     cplex.addRows(lhs, A, rhs)
         %     cplex.addRows(lhs, A, rhs, rowname)
         %
         %   For an example, refer to lpex1.m in examples directory
         %
         %   Parameters:
         %    lhs     double column vector
         %            Lefthand side term for each constraint to be
         %            added to the Cplex problem object.
         %    A       double matrix
         %            The constraint matrix to be added to the Cplex problem
         %            object by rows, it is optional.
         %    rhs     double column vector
         %            Righthand side term for each constraint to be added
         %            to the Cplex problem object, it is a column vector.
         %    rowname char matrix
         %            Names of the new rows, it is optional.
         %
         %   See also Cplex
      end %end of function
      function out = addCols(cpx, obj, A, lb, ub, ctype, colname)
         % Cplex.addCols  Adds columns to a specified Cplex problem object.
         %
         %   Example:
         %    cplex.addCols(obj)
         %    cplex.addCols(obj,A)
         %    cplex.addCols(obj,A,lb)
         %    cplex.addCols(obj,A,lb,ub)
         %    cplex.addCols(obj,A,lb,ub,ctype)
         %    cplex.addCols(obj,A,lb,ub,ctype,name)
         %
         %   For an example, refer to lpex1.m in examples directory
         %
         %   Parameters:
         %    obj     double column vector
         %            Objective function coefficients of the new
         %            variables.
         %    A       double matrix
         %            Constraint matrix to be added to the Cplex problem
         %            object by columns.
         %            Optional.
         %    lb      double column vector
         %            Lower bound on each of the new variables.
         %            Optional.
         %            if lb is not provided,
         %            lb defaults to zeros(length(obj),1)
         %    ub      double column vector
         %            Upper bound on each of the new variables.
         %            Optional.
         %            if ub is not provided ub defaults to
         %            ones(length(obj),1)*Inf
         %    ctype   String
         %            Type of each column in the constraint matrix,
         %            specifying whether a variable is continuous,
         %            integer, binary, semi-continuous, or
         %            semi-integer.
         %            Optional.
         %            if ctype is not provided ctype defaults to
         %            char(ones([1 length(obj)])*('C'))
         %    name    char matrix
         %            Names of the new columns.
         %            Optional.
         %   See also Cplex.
      end %end of function
      function out = delRows(cpx, which)
         % Cplex.delRows  Removes rows with indices listed in 'which' from the
         %   Model.
         %
         %   Example:
         %     cplex.delRows (which)
         %     cplex.delRows ([1 2 5])
         %       deletes the constraints in positions 1, 2 and 5
         %     cplex.delRows ([1:100])
         %       deletes the first 100 constraints
         %
         %   Parameters:
         %     which   double vector
         %          Indices of the constraints (rows) to be
         %          removed.
         %
         %   See also Cplex.
      end %end of function
      function out = delCols(cpx, which)
         % Cplex.delCols  Removes columns with indices listed in
         % 'which' from the Model.
         %
         %   Example:
         %     cplex.delCols (which)
         %     cplex.delCols ([1 2 5])
         %       deletes the variables in positions 1, 2 and 5
         %     cplex.delCols ([1:100])
         %       deletes the first 100 variables
         %    Parameters:
         %     which   double vector
         %          Indices of the constraints (rows) to be
         %          removed.
         %
         %   See also Cplex.
      end %end of function
      function out = addSOSs(cpx, type, ind, wt, varargin)
         % Cplex.addSOSs  Adds information about a special ordered set
         % (SOS) to a problem object of type MILP, MIQP, or MIQCP.
         %  The problem may already contain SOS information.
         %
         % Example
         %   cplex.addSOSs (type, ind, wt, name);
         %   Add single SOS:
         %   cplex.addSOSs ('1', [1 2 3]', [1 2 3]', {'sos1(1)'});
         %   Add set of 2 SOSs:
         %   cplex.addSOSs ...
         %   ('12', {[1 2 3]' [4 5 6]'}, {[2 3 4]' [5 6 7]'}, {'s1' 's2'})
         %   Add SOS structure or structure vector:
         %   cplex.addSOSs(sos);
         %   where sos is structure or structure vector with the fields
         %   type, ind, wt and an optional name.
         %   For an example, refer to mipex3.m in examples directory
         %
         %   Parameters:
         %    type   char '1' or '2' or string
         %         SOS type information for the sets to be added.
         %    ind    double column vector or cell vector
         %         Indices for the sets to be added.
         %    wt     double column vector or cell vector
         %         Weights for the sets to be added.
         %    name   cell vector
         %         Names of the new SOSs. Optional.
         %
         %   See also Cplex.
      end %end of function
      function out = addQCs(cpx, a, Q, sense, rhs, varargin)
         % Cplex.addQCs  Adds a quadratic constraint or quadratic constraint
         % set to a specified Cplex problem object.
         %
         %   Example:
         %     cplex.addQCs(a, Q, sense, rhs, name);
         %     Add single qc
         %     cplex.addQCs([0 0 0]',[1 0 0;0 1 0;0 0 1],'L',1,{'qc1'});
         %     Add qc set(3 qcs)
         %     cplex.addQCs(
         %        [0 0 0;1 0 0;0 1 0]',
         %        {[1 0 0;0 1 0;0 0 1]
         %         [1 0 0;0 1 0;0 0 1]
         %         [1 0 0;0 1 0;0 0 1]},
         %        'LLG',
         %        [1.0 2.0 3.0],
         %        {'qc1' 'qc2' 'qc3'});
         %   Add qc structure or qc structure vector
         %   cplex.addQCs(qc);
         %   qc is structure or structure vector with a, Q, sense,
         %   rhs,and optional name fields.
         %
         %   For an example, refer to qcpex1.m in examples directory.
         %
         %   Parameters:
         %    a     double column vector or matrix
         %         Linear part of the quadratic constraint to be added by
         %         column.
         %    Q     double matrix or cell vector
         %         Quadratic part of the quadratic constraint to be added.
         %    sense   char 'L' or 'G' or string
         %         Sense of the constraint to be added. Note that
         %         quadratic constraints may only be less-than-or-equal-to
         %         greater-than-or-equal-to constraints.
         %    rhs    double or double row vector
         %         Righthand side term for the constraint to be added.
         %    name   cell vector
         %         The name of the constraint to be added. Optional.
         %
         %   See also Cplex.
      end %end of function
      function out = addIndicators(cpx, variable, complemented, a, sense, rhs, varargin)
         % Cplex.addIndicators  Adds indicator constraints to 
         % the specified problem.
         %
         %   Example:
         %     cplex.addIndicators(variable, complemented,
         %                a, sense, rhs, name)
         %   Add one indicator with a name
         %     cplex.addIndicators(5, 1,
         %               [1 2 3 4 5 6]', 'E', 100, {'indc3'});
         %   Add two indicators with names
         %     cplex.addIndicators ([5 6], [0 0],
         %                {[1 1 2 2]' [1 2 3 4 5]'}, 'EE', [10 100],
         %                {'indc3' 'indc4'});
         %   Add indicator structure or structure vector
         %     cplex.addIndicators(indicator);
         %   where indicator is structure or structure vector which has
         %   the fields variable, complemented, a, sense, rhs and
         %   an optional name.
         %
         %   Parameters:
         %    variable double or vector
         %         Binary variable that acts as the indicator for this
         %         constraint.
         %    complemented   double or vector
         %         Boolean value that specifies whether the indicator
         %         variable is complemented.
         %         The linear constraint must be satisfied when the indicator
         %         takes a value of 1 (one) if the indicator is not
         %         complemented, and similarly, the linear constraint must be
         %         satisfied when the indicator takes a value of 0 (zero) if
         %         the indicator is complemented.
         %    a     double column vector or cell vector
         %         Linear portion of the indicator constraint.
         %    sense   char 'L','G' or 'E' or string
         %         Sense of the linear portion of the indicator
         %         constraint. Specify 'L' for <= or 'G' for >= or 'E' for
         %         ==.
         %    rhs    double or double row vector
         %         Righthand side value for the linear portion of the
         %         indicator constraint.
         %    name   optional parameter cell vector
         %         Name of the constraint to be added.
         %         May be empty, in which case the new constraint is assigned a
         %         default name if the indicator constraints already resident
         %         in the Cplex problem object have names; otherwise, no name
         %         is associated with the constraint.
         %
         %   See also Cplex.
      end% end of function
      function out = subsref(cpx, s)
         % Cplex.subsref  Overwrites built-in subsref from MATLAB
         %   See also Cplex.
      function out = subsasgn(cpx,index,val)
         % Cplex.subsasgn  Overwrites built-in subsasgn from MATLAB
         %   See also Cplex.
      end% end of function
   end %end of methods
       function out = loadobj(cpx)
          %   See also Cplex.
       end % end of function
   end %end of methods
      function out = saveobj(cpx)
          %   See also Cplex.
      end % end of function
      function out = defaultNames(cpx, start, num, prefix)
         % The routine defaultNames get the default names used for
         % addCols and addRows.
         %   See also Cplex.
      end % end of function
      function out = isConsistent(cpx)
         % The routine isConsistent check verifies consistency of the
         % Model data and throws exceptions.
         %   See also Cplex.
      end%end of function
      function out = checkStart(cpx)
         % The routine checkStart: checks the Start and MipStart for
         % optimization. There are two kinds of Start, one for continous,
         % another for discrete.  Please refer to the CPLEX User's Manual.
         % Note: the indices here is MipStart from 1 to n.
         %   See also Cplex.
      end%end of function
      function out = isFeasOptConsistent(cpx, preflhs, prefrhs, preflb, prefub)
         % The routine isFeasOptConsistent check the parameters of feasOpt:
         % at least one component of input
         % should be provided.
         %   See also Cplex.
      end %end of function
      function out = getCallbackX(cpx, cbhandle)
         %   See also Cplex.
      end %end of function
      function out = clearSolution(cpx)
         %   See also Cplex.
      end%end of function
      function out = runDisplayFunc(cpx,func_handle,string)
         %   See also Cplex.
      end%end of function
   end%end of methods
end%end of class
