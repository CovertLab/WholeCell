%ReactionProcessTestCase
% Base class for whole cell reaction process test classes. Currently does
% not provide anything beyond the ProcessTestCase class.
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affiliation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 7/13/2010
classdef ReactionProcessTestCase < edu.stanford.covert.cell.sim.ProcessTestCase
    %constructor
    methods
        function this = ReactionProcessTestCase(name)
            this = this@edu.stanford.covert.cell.sim.ProcessTestCase(name);
        end
    end
end