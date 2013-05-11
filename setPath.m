function setPath()
%SETPATH Adds folders to MATLAB path and MATLAB Java path.
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 6/8/2011

if isdeployed
    return;
end

%% MATLAB path

%libraries
setPathHelper('m2html',       'lib/m2html-1.5');                %documentation                   http://www.artefact.tk/software/matlab/m2html/
setPathHelper('propertygrid', 'lib/propertygrid');              %property grid                   http://www.mathworks.com/matlabcentral/fileexchange/28732
setPathHelper('uitabpanel',   'lib/uitabpanel');                %tab panel                       http://www.mathworks.com/matlabcentral/fileexchange/11546
setPathHelper('figstate',     'lib/figstate');                  %min/max figures                 http://www.mathworks.com/support/solutions/en/data/1-3MY8PN/?solution=1-3MY8PN
setPathHelper('spanFigures',  'lib/spanFigures');               %tile figures                    http://www.mathworks.com/matlabcentral/fileexchange/31604
setPathHelper('boundedline',  'lib/boundedline');               %line plot with shaded interval  http://www.mathworks.com/matlabcentral/fileexchange/27485
setPathHelper('tick2text',    'lib/tick2text');                 %tick formatting                 http://www.mathworks.com/matlabcentral/fileexchange/16003
setPathHelper('absolutepath', 'lib/absolutepath');              %full file name                  http://www.mathworks.com/matlabcentral/fileexchange/3857
setPathHelper('runtests',     'lib/matlab_xunit_3.0.1/xunit');  %unit testing                    http://www.mathworks.com/matlabcentral/fileexchange/22846
setPathHelper('multiprod',    'lib/multiprod');                 %tensor multiplication           http://www.mathworks.com/matlabcentral/fileexchange/8773
setPathHelper('glpk',         'lib/glpkmex-2.9');               %linear programming solver       http://glpkmex.sourceforge.net/ (Note win64 version is actually v2.7; maci version is actually v2.8)
setPathHelper('substr',       'lib/util/strutil');              %string utility functions        http://home.online.no/~pjacklam/matlab/software/util/
setPathHelper('iseven',       'lib/util/matutil');              %matrix utility functions        http://home.online.no/~pjacklam/matlab/software/util/
setPathHelper('bullseye',     'lib/bullseye');                  %bullseye polar plot             http://www.mathworks.com/matlabcentral/fileexchange/16458
setPathHelper('swtest',       'lib/swtest');                    %Shapiro-Wilk normality test     http://www.mathworks.com/matlabcentral/fileexchange/13964
setPathHelper('rude',         'lib/rude');                      %run-length encoder/decoder      http://www.mathworks.com/matlabcentral/fileexchange/6436
setPathHelper('Funct_Bezier', 'lib/cubicBezier');               %cubic bezier                    http://www.mathworks.com/matlabcentral/fileexchange/6661
setPathHelper('randsample',   'lib/randsample');                %randsample                      version from R2011b (earlier version has a bug)
setPathHelper('cplexlp',      'lib/cplex-12.2')                 %clpex                           http://www-01.ibm.com/software/websphere/products/optimization/academic-initiative/
setenv('ILOG_LICENSE_FILE', sprintf('%s%s%s%s%s%s%s%s', pwd, filesep, 'lib', filesep, 'cplex-12.2', filesep, 'access.ilm'));
if isempty(which('mxlpsolve'))                                  %lp_solve                        http://sourceforge.net/projects/lpsolve/
    addpath(fullfile(pwd, 'lib/lp_solve_5.5.2.0'));
    switch computer
        case 'PCWIN64', addpath(fullfile(pwd, 'lib\lp_solve_5.5.2.0\bin\win64'));
        case 'PCWIN', addpath(fullfile(pwd, 'lib\lp_solve_5.5.2.0\bin\win32'));
        case {'GLNXA64', 'GLNX86', 'MACI'}, addpath(fullfile(pwd, 'lib/lp_solve_5.5.2.0/bin'));
    end
end

%source code
setPathHelper('edu.stanford.covert.cell.sim.Simulation_Test', 'src_test');  %whole cell test classes
setPathHelper('edu.stanford.covert.cell.sim.Simulation',      'src');       %whole cell classes
addpath(pwd);

%% MATLAB Java path
%source code
setJavaPathHelper('src/+edu/+stanford/+covert/+db/');

%libraries
setJavaPathHelper('lib/json-marshaller/json-0.21.jar');                                 %JSON parser/formatter  http://code.google.com/p/jsonmarshaller/
setJavaPathHelper('lib/mysql-connector-java-5.1.6/mysql-connector-java-5.1.6-bin.jar'); %database connector     http://www.mysql.com/products/connector/
end

%add folder to path, but only if necessary to avoid clearing break points
function setPathHelper(method, pathStr)

%detect if path is already added to MATLAB path, but without calling the
%path method which clear break points
tmp = which(method);
if nargin > 1 && isequal(tmp, [fullfile(pwd, pathStr) pathsep method])
    return;
end

addpath(fullfile(pwd, pathStr));
end

%set java class path, but only if necessary to prevent conflicts with break
%points
function setJavaPathHelper(pathStr)
tmp = javaclasspath;
if pathStr(end) == '/' || pathStr(end) == '\'
    pathStr = pathStr(1:end-1);
end
if any(strcmp(absolutepath(pathStr), tmp))
    return;
end
javaaddpath(absolutepath(pathStr));
end