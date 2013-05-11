%MySQL Database test cases
% Note: the "test" schema and procedures are defined in
% src_test/+edu/+stanford/+covert/+db/test.sql
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affiliation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 2/14/2012
classdef MySQLDatabase_Test < TestCase
    properties
        database
    end
    
    methods
        function this = MySQLDatabase_Test(name)
            this = this@TestCase(name);
        end
    end
    
    methods
        function setUp(this)
            %create database connection
            dbConnectionParameters = config;
            this.database = edu.stanford.covert.db.MySQLDatabase(...
                dbConnectionParameters.hostName, 'test', 'test', 'test');
            this.database.setNullValue(0);
        end
        
        function tearDown(this)
            %close database connection
            this.database.close();
        end
    end
    
    methods
        function testBlobStoringLoading(this)
            %write data to file
            data = char((0:32)');
            fname = tempname;
            fid = fopen(fname,'wb');
            fwrite(fid, data);
            fclose(fid);
            
            %store blob
            this.database.prepareStatement('CALL testBlobIn("{Si}","{F}")', 10001, fname);
            this.database.query();
            delete(fname);
            
            %get last insert id
            this.database.lastInsertID();
            
            %load blob
            this.database.prepareStatement('CALL testBlobOut("{Si}")', 10001);
            result = this.database.query();
            assertEqual(data, char(result.data{1}));
        end
        
        function testReopen(this)
            this.database.reopen();
        end
        
        function testIsClosed(this)
            assertFalse(this.database.isClosed());
        end
        
        function testIsValid(this)
            assertTrue(this.database.isValid());
        end
    end
end