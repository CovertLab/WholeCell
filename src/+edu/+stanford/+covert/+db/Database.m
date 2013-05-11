% Interface for querying a database
%
% Author: Jonathan Karr
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 1/21/2011
classdef Database < handle
    properties (Abstract = true, SetAccess = protected)
        hostName
        schema
        userName
        password
    end

    methods
        function this = Database(varargin)
            switch nargin 
                case 1
                    this.hostName = varargin{1}.hostName;
                    this.schema = varargin{1}.schema;
                    this.userName = varargin{1}.userName;
                    this.password = varargin{1}.password;
                case 4
                    this.hostName = varargin{1};
                    this.schema = varargin{2};
                    this.userName = varargin{3};
                    this.password = varargin{4};
                otherwise
                    throw(MException('Database:error', 'invalid options'));
            end                       
        end
    end

    methods (Abstract = true)
        result = query(this)
    end    
end