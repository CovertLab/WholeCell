%BitMill
% Wrapper for bitmill-bash. Provides several functions:
% - post: submits job to bitmill queue
% - cancel: removes job from bitmill queue
% - getStatus: returns job status(es)
% 
% See also: https://github.com/Numerate/bitmill-bash
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 5/12/2013
classdef BitMill
    properties (Constant = true)
        s3Account = 'prod@numerate.com'
    end
    
    %BitMill methods
    methods (Static = true)
        %Posts job to BitMill queue
        %
        %Input:
        %- type [char]: BitMill job type
        %- parameters [struct]: struct array with two fields -- name and value
        %- inputs [struct]: struct array with two fields -- name and url --
        %  representing the input files
        %- outputs [struct]: struct array with two fields -- name and url --
        %  representing the output files
        %
        %Output:
        %- jobId [char]: BitMill ID of new job
        %- status [double]: 0==>success, otherwise==>failure
        %- errMsg [char]: error message from BitMill
        function [jobId, status, errMsg] = post(varargin)
            import com.numerate.bitmill.BitMill;
            import com.numerate.bitmill.s3cmd;
            
            %% parse inputs
            ip = inputParser();
            
            ip.addParamValue('type', [], @(x) ischar(x));            
            ip.addParamValue('parameters', [], @(x) isstruct(x) && isequal(sort(fieldnames(x)), {'name'; 'value'}));
            ip.addParamValue('inputs', [], @(x) isstruct(x) && isequal(sort(fieldnames(x)), {'name'; 'url'}));
            ip.addParamValue('outputs', [], @(x) isstruct(x) && isequal(sort(fieldnames(x)), {'name'; 'url'}));
            
            ip.parse(varargin{:});
            
            type = ip.Results.type;
            parameters = ip.Results.parameters;
            inputs = ip.Results.inputs;
            outputs = ip.Results.outputs;
            
            validateattributes(type, {'char'}, {'nonempty'});
            
            for i = 1:numel(parameters)
                if isnumeric(parameters(i).value)
                    parameters(i).value = num2str(parameters(i).value);
                elseif ~ischar(parameters(i).value)
                    throw(MException('BitMill:error', 'Invalid parameter value'));
                end
            end
            
            for i = 1:numel(inputs)            
                inputs(i).url = BitMill.s3_to_url(inputs(i).url);
            end
            
            for i = 1:numel(outputs)
                outputs(i).url = BitMill.s3_to_url(outputs(i).url);
            end
            
            %% Run command
            job = struct();
            job.type = type;
            job.parameters = num2cell(parameters(:));
            job.inputs = num2cell(inputs(:));
            job.outputs = num2cell(outputs(:));
            
            cmd = sprintf('post.sh ''%s''', strrep(strrep(...
                savejson('', job, 'ForceRootName', false, 'ArrayIndent', false), ...
                sprintf('\n'), ''), sprintf('\t'), ''));
            [result, status, errMsg] = BitMill.execCmd(cmd);
            jobId = [];
            if status == 0
                tmp = regexp(result, ' id: (?<jobId>.{36,36}) at ', 'names');
                if ~isempty(tmp)
                    jobId = tmp.jobId;
                else
                    status = 1;
                end
            end
        end
        
        %Cancels BitMill job
        %
        %Input:
        %- jobId [string]: BitMill job id
        %
        %Output:
        %- result [char]: former status of canceled job
        %- status [double]: 0==>success, otherwise==>failure
        %- errMsg [char]: error message from BitMill
        function [result, status, errMsg] = cancel(varargin)
            import com.numerate.bitmill.BitMill;
            
            %% parse inputs
            ip = inputParser;
            ip.addRequired('jobId', @(x) ischar(x));
            ip.parse(varargin{:});
            jobId = ip.Results.jobId;
            
            %% Run command
            cmd = sprintf('cancel.sh %s', jobId);
            [result, status, errMsg] = BitMill.execCmd(cmd);
            if status == 0
                if numel(result) > 25
                    status = 0;
                    result = result(25:end-1);
                else
                    status = 1;
                    errMsg = result;
                    result = [];
                end
            end
        end
        
        %Returns BitMill job statuses. If jobId is specified returns status
        %of just the 1 job. Otherwise returns statuses of all jobs.
        %
        %Optional input:
        %- jobId [string]: bitmill job id
        %
        %Output:
        %- result [struct]: struct containing job id, pool, submission time, status, and details
        %- status [double]: 0==>success, otherwise==>failure
        %- errMsg [char]: error message from BitMill
        function [result, status, errMsg] = getStatus(varargin)
            import com.numerate.bitmill.BitMill;
            
            %% parse inputs
            ip = inputParser;
            ip.addOptional('jobId', [], @(x) ischar(x));
            ip.parse(varargin{:});
            jobId = ip.Results.jobId;
            
            %% Run command
            if isempty(jobId)
                cmd = 'jobs.sh';
                [result, status, errMsg] = BitMill.execCmd(cmd);
                if status == 0
                    result = loadjson(result);
                    if isstruct(result) && isfield(result, 'jobs')
                        result = [result.jobs{:}]';
                    else
                        result = repmat(struct('id', [], 'pool', [], 'submitted', [], 'status', [], 'details', []), 0, 1);
                    end
                end
            else
                cmd = sprintf('jobs.sh %s', jobId);
                [result, status, errMsg] = BitMill.execCmd(cmd);
                if status == 0
                    result = loadjson(result);
                end
            end
            
            for i = 1:numel(result)
                result(i).details.parameters = [result(i).details.parameters{:}]';
                result(i).details.inputs = [result(i).details.inputs{:}]';
                result(i).details.outputs = [result(i).details.outputs{:}]';
            end
        end
    end
    
    %helper methods
    methods (Static = true)
        function [result, status, errMsg] = execCmd(cmd)
            if ispc
                cmd = sprintf('bash.exe --login -c "%s"', strrep(cmd, '"', '\"'));
            end
            [status, msg] = system(cmd);
            result = [];
            errMsg = [];
            if status == 0
                result = msg;
            else
                errMsg = msg;
            end
        end
        
        function url = s3_to_url(s3)
            url = regexprep(s3, '^s3:\/\/([^\/]*)\/*', 'https://$1.s3.amazonaws.com/');
        end
    end
end