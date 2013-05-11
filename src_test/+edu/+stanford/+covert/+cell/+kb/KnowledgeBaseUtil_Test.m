classdef KnowledgeBaseUtil_Test < TestCase
    properties
        database
    end
    
    methods
        function this = KnowledgeBaseUtil_Test(name)
            this = this@TestCase(name);
        end
        
        function setUp(this)
            this.database = edu.stanford.covert.db.MySQLDatabase(config);
        end
        
        function tearDown(this)
            this.database.close();
        end
    end
    
    methods
        function testSelectLatestKnowledgeBase(this)
            import edu.stanford.covert.cell.kb.KnowledgeBaseUtil;
            kbWID = KnowledgeBaseUtil.selectLatestKnowledgeBase(this.database);
            validateattributes(kbWID, {'numeric'}, {'positive', 'real', 'integer', 'finite'});
        end
        
        function testSelectLatestSimulation(this)
            import edu.stanford.covert.cell.kb.KnowledgeBaseUtil;
            
            kbWID = KnowledgeBaseUtil.selectLatestKnowledgeBase(this.database);
            simWID = KnowledgeBaseUtil.selectLatestSimulation(kbWID, this.database);
            
            if isempty(simWID)
                simWID2 = this.saveSimulation(kbWID, this.database);
                simWID = KnowledgeBaseUtil.selectLatestSimulation(kbWID, this.database);
                this.deleteSimulation(simWID2, this.database);
                assertEqual(simWID2, simWID);
            end
            
            validateattributes(simWID, {'numeric'}, {'positive', 'real', 'integer', 'finite', 'scalar'});
        end
        
        function testGetSimulationCodeArchive(this)
            import edu.stanford.covert.cell.kb.KnowledgeBaseUtil;
            
            kbWID = KnowledgeBaseUtil.selectLatestKnowledgeBase(this.database);
            simWID = KnowledgeBaseUtil.selectLatestSimulation(kbWID, this.database);
            
            %if no simulation already saved, save small simulation to database
            newSim = 0;
            if isempty(simWID)
                newSim = 1;
                simWID = this.saveSimulation(kbWID, this.database);
            end
            
            %get code archive
            fileName = 'output/runSmallTests/differencesFromRevision.txt';
            if exist(fileName, 'file')
                delete(fileName);
            end
            KnowledgeBaseUtil.getSimulationCodeArchive(simWID, fileName, this.database);
            assertEqual(2, exist(fileName, 'file'));
            
            %delete the added small simulation
            if newSim
                this.deleteSimulation(simWID, this.database);
            end
        end
        
        function testGetContact(this)
            import edu.stanford.covert.cell.kb.KnowledgeBaseUtil;
            
            contact = struct('email', 'jkarr@stanford.edu');
            
            kbWID = KnowledgeBaseUtil.selectLatestKnowledgeBase(this.database);
            wid = KnowledgeBaseUtil.getContact(contact, kbWID, this.database);
            validateattributes(wid, {'numeric'}, {'positive', 'real', 'integer', 'finite'});
        end
        
        function testNewContact(this)
            import edu.stanford.covert.cell.kb.KnowledgeBaseUtil;
            
            contact = struct(...
                'firstName', 'Jonathan', ...
                'lastName', 'Karr', ...
                'affiliation', 'Covert Lab, Department of Bioengineering, Stanford University', ...
                'email', 'jkarr@stanford.edu');
            
            kbWID = KnowledgeBaseUtil.selectLatestKnowledgeBase(this.database);
            wid = KnowledgeBaseUtil.newContact(contact, kbWID, this.database);
            validateattributes(wid, {'numeric'}, {'positive', 'real', 'integer', 'finite'});
        end
    end
    
    %helper methods
    methods (Static)
        function simWID = saveSimulation(kbWID, database)
            import edu.stanford.covert.cell.sim.SimulationFixture;
            import edu.stanford.covert.cell.sim.util.ConditionSet;
            import edu.stanford.covert.cell.sim.util.DatabaseLogger;
            
            sim = SimulationFixture.load([], true);
            sim.applyOptions('lengthSec', 10);
            
            metadata = ConditionSet.parseConditionSet(sim, 'src_test/+edu/+stanford/+covert/+cell/+sim/+util/fixtures/runSimulations.xml').metadata;
            metadata.knowledgeBaseWID = kbWID;
            metadata.revision = edu.stanford.covert.util.revision;
            
            downsampleStepSec = 2;
            logger = DatabaseLogger(database, downsampleStepSec, metadata);
            logger.addMetadata(metadata);
            sim.run(logger);
            
            simWID = logger.simulationWID;
        end
        
        function deleteSimulation(simWID, database)
            database.prepareStatement('CALL delete_simulation("{Si}")', simWID);
            database.query();
        end
    end
end