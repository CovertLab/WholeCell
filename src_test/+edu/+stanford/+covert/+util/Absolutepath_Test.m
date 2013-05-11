classdef Absolutepath_Test < TestCase
    methods
        function this = Absolutepath_Test(name)
            this = this@TestCase(name);
        end
        
        function testAbsolutepath(~)
            tmp = pwd;
            if ispc; tmp(1) = upper(tmp(1)); end;
            assertEqual(tmp(1:find(tmp==filesep,1,'last')-1),  absolutepath('..'));
            assertEqual(tmp,  absolutepath('.'));
            assertEqual([tmp filesep 'runSmallTests.m'],  absolutepath('./runSmallTests.m'));
            assertEqual([tmp filesep 'runSmallTests.m'],  absolutepath('../simulation/./runSmallTests.m'));
            assertEqual([tmp filesep 'runSmallTests.m'],  absolutepath('../simulation/../simulation/./runSmallTests.m'));
            assertEqual([tmp filesep 'output' filesep 'runSmallTests'],  absolutepath('../simulation/../simulation/./output/runSmallTests'));
            assertEqual([tmp filesep 'output' filesep 'runSmallTests'],  absolutepath('..//simulation/../simulation/./output/runSmallTests'));
            tic; absolutepath('../simulation/../simulation/./output/runSmallTests'); toc
        end
    end
end