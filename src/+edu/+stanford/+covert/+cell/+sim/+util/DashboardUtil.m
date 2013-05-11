classdef DashboardUtil
    properties (Constant)
        dashboardURL = 'http://covertlab.stanford.edu/projects/wholecell/browser/bin/index.php';
    end
    
    methods (Static = true)
        %opens the dashboard in browser
        function openDashboard(knowledgeBaseWID, simulationWID)
            url = edu.stanford.covert.cell.sim.util.DashboardUtil.dashboardURL;
            if exist('knowledgeBaseWID','var')
                url = [url '?knowledgeBaseWID=' num2str(knowledgeBaseWID)];
                if exist('simulationWID','var')
                    url = [url '&simulationWID=' num2str(simulationWID)];
                end
            end
            web(url, '-browser');
        end
    end
end