%s3cmd
% Wrapper for s3cmd. Provides several functions:
% - makeBucket: creates bucket
% - listBuckets: lists available S3 buckets
% - put: saves file to S3 bucket
% - get: downloads file from S3 bucket
% - delete: deletes file
% - listContents: lists contents of S3 bucket
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 5/12/2013
classdef s3cmd
    methods (Static = true)
        function [status, errMsg] = makeBucket(bucketUrl)
            import com.numerate.bitmill.s3cmd;
            [~, status, errMsg] = s3cmd.execCmd(sprintf('s3cmd mb %s', bucketUrl));
        end
        
        function [buckets, status, errMsg] = listBuckets()
            import com.numerate.bitmill.s3cmd;
            [result, status, errMsg] = s3cmd.execCmd('s3cmd ls');
            if status == 0
                result = result(:, 1:end-1);
                
                nBuckets = size(result, 1);
                buckets = repmat(struct('url', [], 'timestamp', []), nBuckets, 1);
                for i = 1:nBuckets
                    buckets(i).url = result(19:find(result ~= ' ', 1, 'last'));
                    buckets(i).timestamp = datenum(result(i, 1:16));
                end
            else
                buckets = [];
            end
        end
        
        function [status, errMsg] = put(localFilePath, remoteFileUrl)
            import com.numerate.bitmill.s3cmd;
            cmd = sprintf('s3cmd put %s %s', s3cmd.escapeLocalPath(localFilePath), remoteFileUrl);
            [~, status, errMsg] = s3cmd.execCmd(cmd);
        end
        
        function [status, errMsg] = get(remoteFileUrl, localFilePath, force)
            import com.numerate.bitmill.s3cmd;
            
            forceOpt = '';
            if nargin >= 3 && force
                forceOpt = '--force';
            end
            
            cmd = sprintf('s3cmd get %s %s %s', forceOpt, remoteFileUrl, s3cmd.escapeLocalPath(localFilePath));
            [~, status, errMsg] = s3cmd.execCmd(cmd);
        end
        
        function [status, errMsg] = delete(fileUrl)
            import com.numerate.bitmill.s3cmd;
            [~, status, errMsg] = s3cmd.execCmd(sprintf('s3cmd del %s', fileUrl));
        end
        
        function [contents, status, errMsg] = listContents(bucketUrl)
            import com.numerate.bitmill.s3cmd;
            [result, status, errMsg] = s3cmd.execCmd(sprintf('s3cmd ls %s', bucketUrl));
            if status == 0
                result = strsplit(sprintf('\n'), result);
                
                nFiles = numel(result) - 1;
                contents = repmat(struct('url', [], 'timestamp', [], 'size', []), nFiles, 1);
                for i = 1:nFiles
                    tmp = strsplit(' ', result{i}(18:end));
                    tmp = tmp(~cellfun(@isempty, tmp));
                    contents(i).url = tmp{2};
                    contents(i).size = str2double(tmp{1});
                    contents(i).timestamp = datenum(result{i}(1:16));
                end
            else
                contents = [];
            end
        end
        
        function [status, errMsg] = grantAclRead(remoteFile, account)
            import com.numerate.bitmill.s3cmd;
            cmd = sprintf('s3cmd setacl --acl-grant=read:%s %s', account, remoteFile);
            [~, status, errMsg] = s3cmd.execCmd(cmd);
        end
        
        function [status, errMsg] = grantAclReadAcp(remoteFile, account)
            import com.numerate.bitmill.s3cmd;
            [~, status, errMsg] = s3cmd.execCmd(sprintf('s3cmd setacl --acl-grant=read_acp:%s %s', account, remoteFile));
        end
        
        function [status, errMsg] = grantAclWrite(remoteFile, account)
            import com.numerate.bitmill.s3cmd;
            [~, status, errMsg] = s3cmd.execCmd(sprintf('s3cmd setacl --acl-grant=write:%s %s', account, remoteFile));
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
        
        %escape for cygwin
        function path = escapeLocalPath(path)
            if ispc
                if path(2) ~= ':'
                    path = fullfile(pwd, path);
                end
                
                path(path == '\') = '/';
                
                path = sprintf('/cygdrive/%s/%s', lower(path(1)), path(4:end));
            end
        end
    end
end