%Animation
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affiliation: Covert Lab, Department of Bioengineering, Stanford University
% Last Updated: 7/20/2011
classdef Animation < handle
    properties (Abstract = true, SetAccess = protected)
        title
        description
        author
    end
    
    properties
        year
        copyright
                
        frameFilePattern = 'frame_%d'
        frameFormat = 'png'
        frameFormatOptions
        subtitleFileName = 'subtitles.srt'
        movieFileName
        movieDPI = 250;
        movieWidth
        movieHeight
        movieOptions
        
        frameRenderer
        movieRenderer
        codec = struct(...
            'ffmpeg', 'libx264', ...
            'mencoder', 'mpeg4')
        subtitleCodec = 'srt'
        frameRate = 12
        quality = 4
        bitRate = 1200
        loopCount = 0
        
        verbosity
    end
    
    properties
        tmpDirectory
    end
    
    properties (SetAccess = protected)
        simGroupId
        simId
        downsampleStepSizeSec
        
        simulation
        times
        states
    end
    
    properties (SetAccess = protected)
        subtitleFid
    end
    
    methods
        function this = Animation(simGroupId, simId, movieFileName, downsampleStepSizeSec, frameRenderer, movieRenderer, verbosity)
            import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;
            
            %% options
            if nargin < 1 || isempty(simGroupId)
                simGroupId = SimulationDiskUtil.getLatestSimulationGroup();
            end
            if nargin < 2 || isempty(simId)
                simId = 1;
            end
            if nargin < 3 || isempty(movieFileName)
                movieFileName = 'video.avi';
            end
            if nargin < 4 || isempty(downsampleStepSizeSec)
                downsampleStepSizeSec = 100;
            end
            if nargin < 5 || isempty(frameRenderer)
                frameRenderer = 'inkscape';
            end
            if nargin < 6 || isempty(movieRenderer)
                movieRenderer = 'ffmpeg';
            end
            if nargin < 7 || isempty(verbosity)
                verbosity = 0;
            end
            
            this.simGroupId = simGroupId;
            this.simId = simId;
            this.movieFileName = movieFileName;
            this.downsampleStepSizeSec = downsampleStepSizeSec;
            this.frameRenderer = frameRenderer;
            this.movieRenderer = movieRenderer;
            this.verbosity = verbosity;
            
            this.year = datestr(now, 'YYYY');
            this.copyright = ['Copyright ' this.author ' ' this.year];
            this.tmpDirectory = tempname;
            
            this.movieWidth = 600;
            this.movieHeight = 900;
            this.frameFormatOptions.png = '';
            this.movieOptions.ffmpeg = '';
            if ispc
                this.frameFormatOptions.png = '-define png:bit-depth=8 -define png:color-type=2 -define png:compression-level=0 -define png:compression-strategy=0';
            else
                this.movieOptions.ffmpeg = '-vpre medium';
            end
            
            %% load simulation and data
            [~, ~, this.simulation] = SimulationDiskUtil.getSimulation([this.simGroupId filesep num2str(this.simId(1))]);
            [this.times, this.states] = this.loadSimulationData();
        end
        
        function run(this, frameSkip, drawFrames, renderFrames, renderMovie, cleanup)
            if nargin < 2 || isempty(frameSkip)
                frameSkip = [1 1];
            end
            
            nFrames = floor(numel(this.times) / this.downsampleStepSizeSec);
            
            %make temporary directory
            if (nargin < 6 || cleanup) && ~exist(this.tmpDirectory, 'dir')
                mkdir(this.tmpDirectory);
            end
            
            %draw frames
            if nargin < 3 || drawFrames
                for frameIdx = frameSkip(1):frameSkip(2):nFrames
                    %calculate time index
                    timeIdx = (frameIdx - 1) * this.downsampleStepSizeSec + 1;
                    
                    %build svg
                    fileName = sprintf(['%s%s' this.frameFilePattern '.svg'], this.tmpDirectory, filesep, frameIdx);
					svg = this.drawFrame(timeIdx);
                    fid = fopen(fileName, 'w');
                    fwrite(fid, svg);
                    fclose(fid);
                end
            end
            
            %convert to png
            if nargin < 4 || renderFrames
                this.renderFrames(frameSkip);
            end
            
            %render movie
            if nargin < 5 || renderMovie
                % subtitles
                this.initializeSubtitles();
                for frameIdx = frameSkip(1):frameSkip(2):nFrames
                    %calculate time index
                    timeIdx = (frameIdx - 1) * this.downsampleStepSizeSec + 1;
                    
                    %append
                    this.appendSubtitle(frameIdx, timeIdx);
                end
                this.finalizeSubtitles();
                
                %movie
                this.renderMovie();
            end
            
            %clean up temporary directory
            if nargin < 6 || cleanup
                rmdir(this.tmpDirectory, 's');
            end
        end
        
        function renderFrames(this, frameSkip)
            switch this.frameRenderer
                case 'librsvg'
                    this.renderFrames_librsvg(frameSkip);
                case 'graphicsmagick'
                    this.renderFrames_graphicsmagick(frameSkip);
                case 'imagemagick'
                    this.renderFrames_imagemagick(frameSkip);
                case 'inkscape'
                    this.renderFrames_inkscape(frameSkip);
                case 'batik'
                    this.renderFrames_batik(frameSkip);
                otherwise
                    throw(MException('Animation:error', 'Undefined renderer %s', this.frameRenderer))
            end
        end
        
        function renderFrames_librsvg(this, frameSkip)
            nFrames = floor(numel(this.times) / this.downsampleStepSizeSec);
            renderFileName = [this.tmpDirectory, filesep 'renderFrames.sh'];
            fid = fopen(renderFileName, 'w');
            fprintf(fid, '#!/bin/sh\n');
            for frameIdx = frameSkip(1):frameSkip(2):nFrames
                if ~exist(sprintf('%s%s%s.svg', this.tmpDirectory, filesep, strrep(this.frameFilePattern, '%d', num2str(frameIdx))), 'file')
                    continue;
                end
                fprintf(fid, [...
                    'rsvg-convert '...
                    '-f png '...
                    '-w %d -h %d '...
                    '-o "%s%s%s.png" "%s%s%s.svg"\n'], ...
                    this.movieWidth, this.movieHeight, ...
                    this.tmpDirectory, filesep, strrep(this.frameFilePattern, '%d', num2str(frameIdx)), ...
                    this.tmpDirectory, filesep, strrep(this.frameFilePattern, '%d', num2str(frameIdx)));
            end
            fclose(fid);
            [status, result] = system(['chmod 775 ' renderFileName ' && ' renderFileName]);
            if status ~= 0
                throw(MException('Animation:error', 'Failed to convert svg to png: %s', result));
            end
        end
        
        function renderFrames_graphicsmagick(this, ~)
            cmd = sprintf([...
                'gm mogrify '...
                '-format png %s '...
                '-resize %dx%d '...
                '"%s%s%s.svg"'], ...
                this.frameFormatOptions.png, this.movieWidth, this.movieHeight, this.tmpDirectory, filesep, strrep(this.frameFilePattern, '%d', '*'));
            [status, result] = system(cmd);
            if status ~= 0
                throw(MException('Animation:error', 'Failed to convert svg to png: %s', result));
            end
        end
        
        function renderFrames_imagemagick(this, ~)
            cmd = sprintf([...
                'mogrify '...
                '-format png %s '...
                '-resize %dx%d '...
                '"%s%s%s.svg"'], ...
                this.frameFormatOptions.png, this.movieWidth, this.movieHeight, this.tmpDirectory, filesep, strrep(this.frameFilePattern, '%d', '*'));
            [status, result] = system(cmd);
            if status ~= 0
                throw(MException('Animation:error', 'Failed to convert svg to png: %s', result));
            end
        end
        
        function renderFrames_inkscape(this, frameSkip)
            nFrames = floor(numel(this.times) / this.downsampleStepSizeSec);
            for frameIdx = frameSkip(1):frameSkip(2):nFrames
                if ~exist(sprintf('%s%s%s.svg', this.tmpDirectory, filesep, strrep(this.frameFilePattern, '%d', num2str(frameIdx))), 'file')
                    continue;
                end
                
                %print status
                if this.verbosity > 0
                    fprintf('Drawing frame %d of %d frames\n', ...
                        frameIdx, nFrames);
                end
                [status, result] = system(sprintf('inkscape "%s" --export-png="%s" --export-area-page -w %d -h %d', ...
                    sprintf(['%s%s' this.frameFilePattern '.svg'], this.tmpDirectory, filesep, frameIdx), ...
                    sprintf(['%s%s' this.frameFilePattern '.png'], this.tmpDirectory, filesep, frameIdx), ...
                    this.movieWidth, this.movieHeight));
                if status ~= 0
                    throw(MException('Animation:error', 'Failed to convert svg to png: %s', result));
                end
            end
        end
        
        function renderFrames_batik(this, ~)
            [status, result] = system(sprintf('java -jar batik-rasterizer.jar -w %d -h %d -m image/png -d "%s" "%s%s%s.svg"', ...
                this.movieWidth, this.movieHeight, this.tmpDirectory, this.tmpDirectory, filesep, strrep(this.frameFilePattern, '%d', '*')));
            if status ~= 0
                throw(MException('Animation:error', 'Failed to convert svg to png: %s', result));
            end
        end
        
        function initializeSubtitles(this)
            this.subtitleFid = fopen([this.tmpDirectory filesep this.subtitleFileName], 'w');
        end
        
        function appendSubtitle(this, frameIdx, timeIdx)
            fprintf(this.subtitleFid, '%d\n', frameIdx);
            fprintf(this.subtitleFid, '%s --> %s\n', ...
                datestr((frameIdx-1)/this.frameRate/(24*60*60), 'HH:MM:SS,FFF'), ...
                datestr((frameIdx-1)/this.frameRate/(24*60*60), 'HH:MM:SS,FFF'));
            fprintf(this.subtitleFid, '%s\n\n', datestr(this.times(timeIdx) / 3600 / 24, 'HH:MM:SS'));
        end
        
        function finalizeSubtitles(this)
            fclose(this.subtitleFid);
        end
        
        function renderMovie(this)
            switch this.movieRenderer
                case 'ffmpeg'
                    this.renderMovie_ffmpeg();
                case 'mencoder'
                    this.renderMovie_mencoder();
                otherwise
                    throw(MException('Animation:error', 'Undefined renderer %s', this.movieRenderer))
            end
        end
        
        function renderMovie_ffmpeg(this)
            if exist(this.movieFileName, 'file')
                delete(this.movieFileName);
            end
            cmd = sprintf([...
                'ffmpeg -f image2 -i "%s/%s.png" ' ...
                '-vcodec %s -r %f %s -qscale %d -qmin 1 -qmax 35 -loop_output %d -threads 0 -s %dx%d ' ...
                '-an ' ...
                '-timestamp now -metadata title="%s" -metadata description="%s" -metadata author="%s" -metadata year="%s"  -metadata copyright="%s" ' ...
                '-y "%s"'
                ], ...
                this.tmpDirectory, this.frameFilePattern, ...
                this.codec.ffmpeg, this.frameRate, this.movieOptions.ffmpeg, this.quality, this.loopCount, this.movieWidth, this.movieHeight, ...
                this.title, this.description, this.author, this.year, this.copyright, ...
                this.movieFileName);
            [status, result] = system(cmd);
            if status ~= 0
                throw(MException('Simulation:error', 'Failed to render movie: %s', result));
            end
        end
        
        function renderMovie_mencoder(this)
            movieFileName = this.movieFileName; %#ok<*PROP>
            if ~((isunix && movieFileName(1) == '/') || (ispc && movieFileName(2) == ':'))
                movieFileName = [pwd filesep movieFileName];
            end
            
            cmd = sprintf([
                'cd %s && ' ...
                'mencoder "mf://%s.png" -mf type=png:fps=%f:w=%d:h=%d -sub "%s" -o "%s" ' ...
                '-ovc lavc -lavcopts vcodec=%s:vbitrate=%d:vpass=1 ' ...
                '-subcp enca:en:iso-8859-7 ' ...
                '-info name="%s":subject="%s":artist="%s":copyright="%s" ' ...
                ], ...
                this.tmpDirectory, this.frameFilePattern, this.frameRate, this.movieWidth, this.movieHeight, ...
                this.subtitleFileName, movieFileName, ...
                this.codec.mencoder, this.bitRate, ...
                this.title, this.description, this.author, this.copyright);
            [status, result] = system(cmd);
            if status ~= 0
                throw(MException('Simulation:error', 'Failed to render movie: %s', result));
            end
        end
    end
    
    methods (Abstract = true)
        [times, states] = loadSimulationData(this)
        
        svg = drawFrame(this, timeIdx)
    end
end