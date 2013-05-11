function setWarnings()
%SETWARNINGS Turns on and off (mostly off) specific warnings.
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 9/21/2010

warning('off', 'MATLAB:dispatcher:nameConflict');
warning('off', 'MATLAB:class:cannotUpdateClass');
warning('off', 'MATLAB:class:cannotUpdateClass:Changed');
warning('off', 'MATLAB:ClassInstanceExists');
warning('off', 'MATLAB:xlsread:Mode')
warning('off', 'MATLAB:xlsread:RangeIncompatible');
warning('off', 'MATLAB:javaclasspath:duplicateEntry');
warning('off', 'MATLAB:xlswrite:AddSheet');
warning('off', 'MATLAB:structOnObject');
warning('off', 'MATLAB:Print:CustomResizeFcnInPrint');
warning('off', 'MATLAB:save:sizeTooBigForMATFile');
warning('off', 'MATLAB:load:variableNotFound');
warning('off', 'WholeCell:warning:network');
warning('off', 'WholeCell:warning:initializeNegativeCounts');
warning('off', 'MATLAB:Axes:NegativeDataInLogAxis');
warning('off', 'swtest:largeSampleSize');
warning('off', 'stats:gamfit:ZerosInData');
warning('off', 'MATLAB:RandStream:loadobj:LoadError');