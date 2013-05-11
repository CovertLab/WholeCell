%Simulation logger interface:
%- initialize
%- append
%- finalize
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 1/10/2011
classdef Logger < handle
    methods (Abstract)
        this = initialize(this, sim)
        this = append(this, sim)
        this = finalize(this, sim)
    end    
end