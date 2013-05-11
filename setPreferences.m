function setPreferences()
%SETPREFERENCES Sets MATLAB preferences
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Author: Derek Macklin, macklin@stanford.edu
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 9/6/2011

%prevent MATLAB from using OpenGL to avoid segmentation faults
set(0, 'DefaultFigureRenderer', 'painters')
opengl('neverselect')

%set character encoding to UTF-8
if ~isdeployed
    slCharacterEncoding('UTF-8');
end

%default font
set(0, 'defaultAxesFontName', 'Arial')
set(0, 'defaultTextFontName', 'Arial')