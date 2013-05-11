%CellGeometryAnimation
%
% Author: Derek Macklin, macklin@stanford.edu
%         (Based on code developed by Jonathan Karr, jkarr@stanford.edu)
% Affiliation: Covert Lab, Department of Bioengineering, Stanford University
% Last Updated: 7/22/2011
classdef CellGeometryAnimation < handle
    properties (SetAccess = protected)
        title = 'Whole Cell DNA Animation'
        description = 'Animation of cell geometry dynamics in a single Mycoplasma genitalium whole-cell simulation.'
        author = 'Derek Macklin, Covert Lab, Department of Bioengineering, Stanford University'
        year
        copyright
        
        tmpDirectory
        frameFilePattern = 'frame_%d'
        subtitleFileName = 'subtitles.srt'
        movieFileName
        movieDPI = 250;
        movieWidth
        movieHeight
        
        frameRenderer
        movieRenderer
        codec = struct(...
            'ffmpeg', 'libx264', ...
            'mencoder', 'mpeg4')
        subtitleCodec = 'srt'
        frameRate = 30
        bitRate = 1200
        loopCount = 0
        
        verbosity
    end
    
    properties (SetAccess = protected)
        simGroupId
        simId
        downsampleStepSizeSec
        
        simulation
        states
        log
        feature
    end
    
    properties (SetAccess = protected)
        subtitleFid
    end
    
    methods
        function this = CellGeometryAnimation(simGroupId, simId, movieFileName, downsampleStepSizeSec, frameRenderer, movieRenderer, verbosity)
            if nargin < 1 || isempty(simGroupId)
                simGroupId = '2011_07_15_20_08_43';
            end
            if nargin < 1 || isempty(simId)
                simId = 1;
            end
            if nargin < 3 || isempty(movieFileName)
                movieFileName = 'video.avi';
            end
            if nargin < 4 || isempty(downsampleStepSizeSec)
                downsampleStepSizeSec = 100;
            end
            if nargin < 5 || isempty(frameRenderer)
                frameRenderer = 'imagemagick';
            end
            if nargin < 6 || isempty(movieRenderer)
                movieRenderer = 'mencoder';
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
            
            this.movieWidth = ceil(600/90 * this.movieDPI);
            this.movieHeight = ceil(300/90 * this.movieDPI);
        end
        
        function run(this)
            %make temporary directory
            mkdir(this.tmpDirectory);
            
            %load simulation data
            this.loadSimulationData();
            
            %initialize subtitles
            this.initializeSubtitles();
            
            %draw frames and build subtitles
            nFrames = floor(numel(this.states.Time.values) / this.downsampleStepSizeSec);
            for frameIdx = 1:nFrames
                %calculate time index
                timeIdx = (frameIdx - 1) * this.downsampleStepSizeSec + 1;
                
                %build svg
                this.drawFrame(frameIdx, timeIdx);
                
                %subtitle
                this.appendSubtitle(frameIdx, timeIdx);
            end
            
            %finalize subtitles
            this.finalizeSubtitles();
            
            %convert to png
            this.renderFrames();
            
            %render movie
            this.renderMovie();
            
            %clean up temporary directory
            rmdir(this.tmpDirectory, 's');
        end
        
        function loadSimulationData(this)
            import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;
            import edu.stanford.covert.cell.sim.util.DiskLogger;
            import edu.stanford.covert.cell.sim.state.CellGeometry;
            [simDir,~,this.simulation]=SimulationDiskUtil.getSimulation('2011_07_19_18_35_32/1');
            stateNames={'Geometry' 'width'
                        'Geometry' 'pinchedDiameter'
                        'Chromosome' 'polymerizedRegions'
                        'Chromosome' 'monomerBoundSites'
                        'Chromosome' 'complexBoundSites'
                        'Time' 'values'};
            this.states = DiskLogger.load(simDir, stateNames, [], [], [], 'extract');
            this.log = load([simDir filesep 'summary.mat']);
            
            cm = this.simulation.state('Mass');
            geom = this.simulation.state('Geometry');
            density = geom.density;
            mass = (this.log.mass(2:end) * cm.cellInitialDryWeight)/(1-cm.fractionWetWeight);
            width = squeeze(this.states.Geometry.width);
            pinchedDiameter = squeeze(this.states.Geometry.pinchedDiameter);
            t = squeeze(this.states.Time.values);
            septumLength = (width - pinchedDiameter)/2;
            
            cylindricalLength = zeros(size(mass));
            maxTotalLengths = zeros(size(mass));
            
            for i=1:length(cylindricalLength)
                if(pinchedDiameter(i) > 0)
                    [cylindricalLength(i), ~, maxTotalLengths(i)] = CellGeometry.calculateGeometry(...
                        width(i), ...
                        pinchedDiameter(i), ...
                        mass(i)/density);
                end
            end
            
            this.feature = struct;
            this.feature.cylindricalLength = cylindricalLength;
            this.feature.septumLength = septumLength;
            this.feature.maxTotalLengths = maxTotalLengths;
            this.feature.t = t;

        end
        
        function drawFrame(this, frameIdx, timeIdx)
            import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;
            import edu.stanford.covert.cell.sim.util.DiskLogger;
            import edu.stanford.covert.cell.sim.state.CellGeometry;
            fileName = sprintf(['%s%s' this.frameFilePattern '.svg'], this.tmpDirectory, filesep, frameIdx);
            chr = this.simulation.state('Chromosome');

            cylindricalLength = this.feature.cylindricalLength;
            septumLength = this.feature.septumLength;
            maxTotalLengths = this.feature.maxTotalLengths;
            t = this.feature.t;
            width = this.states.Geometry.width;
            
            scale = 1e-9;
            
            idx = find(t==timeIdx);
            maxTotalLength = max(maxTotalLengths);
            
            svg_width = maxTotalLength / scale + 2;
            svg_height = max(width) / scale + 2;
            
            %open svg
            fid = fopen(fileName, 'w');
            fprintf(fid, '<?xml version="1.0" standalone="no"?>\n');
            fprintf(fid, '<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN" "http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">\n');
            fprintf(fid, '<svg width="%d" height="%d" version="1.1" xmlns="http://www.w3.org/2000/svg" viewBox="%d %d %d %d">\n', svg_width, svg_height, 0, 0, svg_width, svg_height);
            fprintf(fid, '<defs>\n');
            fprintf(fid, '<linearGradient id="cellBg" x1="0%%" y1="0%%" x2="0%%" y2="100%%">\n');
            fprintf(fid, '<stop style="stop-color:#97b0fb;stop-opacity:0.05" offset="0"/>\n');
            fprintf(fid, '<stop style="stop-color:#89b8d6;stop-opacity:1" offset="1"/>\n');
            fprintf(fid, '</linearGradient>');
            fprintf(fid, '</defs>\n');

            % White background
            fprintf(fid, '<rect x="%d" y="%d" width="%d" height="%d" style="fill:white; stroke:none"/>\n', ...
                0, 0, svg_width, svg_height);
           
            % Draw the cell
            fprintf(fid, this.drawCell(scale,width(idx),cylindricalLength(idx),septumLength(idx),maxTotalLength,max(width)));
            
            % Draw the chromosomes
            [posStrnds, lens] = find(this.states.Chromosome.polymerizedRegions);
            lens = lens(posStrnds(:, 3) == timeIdx, :);
            posStrnds = posStrnds(posStrnds(:, 3) == timeIdx, 1:2);
            centerX = [ (maxTotalLength / 2 - septumLength(idx) - cylindricalLength(idx) / 2), (maxTotalLength / 2 + septumLength(idx) + cylindricalLength(idx) / 2) ];
            centerY = [ (width(idx) / 2), (width(idx) / 2) ];
            fprintf(fid, this.drawArc(posStrnds, lens, [width(idx)*0.5 width(idx)*0.55]/(2*scale), 'black', 2, this.log.replisome(1:2,idx), centerX/scale, centerY/scale));
            
            %bound proteins (RNA Pol)
            [cpxPosStrnds, cpxIdxs] = find(this.states.Chromosome.complexBoundSites == 200);
            posStrnds = [...
                cpxPosStrnds(cpxPosStrnds(:, 3) == timeIdx, 1:2)
                ];
            lens = [...
                chr.complexDNAFootprints(cpxIdxs(cpxPosStrnds(:, 3) == timeIdx, :))
                ];
            fprintf(fid, this.drawCircles(posStrnds, lens, [width(idx)*0.45 width(idx)*0.60]/(2*scale), 'green', 1, this.log.replisome(1:2,idx), centerX/scale, centerY/scale, 5, 'none')); 

            clear cpxPosStrnds;
            clear cpxIdxs;
            %bound proteins (Helicase)
            [cpxPosStrnds, cpxIdxs] = find(this.states.Chromosome.complexBoundSites == 50);
            if numel(cpxIdxs) > 0
                posStrnds = [...
                    cpxPosStrnds(cpxPosStrnds(:, 3) == timeIdx, 1:2)
                    ];
                lens = [...
                    chr.complexDNAFootprints(cpxIdxs(cpxPosStrnds(:, 3) == timeIdx, :))
                    ];
                if ~isempty(posStrnds)
                    fprintf(fid, this.drawCircles(posStrnds, lens, [width(idx)*0.45 width(idx)*0.60]/(2*scale), 'gold', 1, this.log.replisome(1:2,idx), centerX/scale, centerY/scale, 3, 'gold'));
                end
            end

            clear cpxPosStrnds;
            clear cpxIdxs;
             %bound proteins (DnaA-ADP 1mer)
            [cpxPosStrnds, cpxIdxs] = find(this.states.Chromosome.complexBoundSites == 181);
            if numel(cpxIdxs) > 0
                posStrnds = [...
                    cpxPosStrnds(cpxPosStrnds(:, 3) == timeIdx, 1:2)
                    ];
                lens = [...
                    chr.complexDNAFootprints(cpxIdxs(cpxPosStrnds(:, 3) == timeIdx, :))
                    ];
                if ~isempty(posStrnds)
                    fprintf(fid, this.drawCircles(posStrnds, lens, [width(idx)*0.525 width(idx)*0.525]/(2*scale), 'lightskyblue', 1, this.log.replisome(1:2,idx), centerX/scale, centerY/scale, 2, 'lightskyblue'));
                end
            end
            
            %close svg
            fprintf(fid, '</svg>\n');
            fclose(fid);
        end
        
        function svg = drawCell(this, scale, width, cylindricalLength, septumLength, maxTotalLength, maxWidth)
            width = width / scale;
            cylindricalLength = cylindricalLength / scale;
            septumLength = septumLength / scale;
            maxTotalLength = maxTotalLength / scale;
            maxWidth = maxWidth / scale;
            
            svg = ['<path d="\n', ...
                   sprintf('M%0.4f,1 a%0.4f,%0.4f 0 1 0 0,%0.4f\n', maxTotalLength/2-septumLength-cylindricalLength/2+1, width/2, width/2, width), ...
                   sprintf('l%0.4f,0\n', cylindricalLength/2), ...
                   sprintf('l%0.4f,-%0.4f\n', septumLength, septumLength), ...
                   sprintf('m0,0 l%0.4f,%0.4f\n', septumLength, septumLength), ...
                   sprintf('l%0.4f,0\n', cylindricalLength/2), ...
                   sprintf('a%0.4f,%0.4f 0 1 0 0,-%0.4f\n', width/2,width/2, width), ...
                   sprintf('l-%0.4f,0\n', cylindricalLength/2), ...
                   sprintf('l-%0.4f,%0.4f\n', septumLength, septumLength), ...
                   sprintf('l-%0.4f,-%0.4f\n', septumLength, septumLength), ...
                   sprintf('l-%0.4f,0"\n', cylindricalLength/2), ...
                   'style="stroke:rgb(0,0,0);stroke-width:2;stroke-linecap:round;stroke-linejoin:miter;fill:url(#cellBg)"/>\n'];
        end
        
        function svg = drawArc(this, posStrnds, lens, r0, strokeColor, strokeWidth, helicase, centerX, centerY)
            chr = this.simulation.state('Chromosome');
            strokeColor_set = strokeColor;
            
            svg = [];
            for j = 1:size(posStrnds, 1)
                theta = posStrnds(j, 1) / chr.sequenceLen * 2 * pi - pi/2;
                dTheta = lens(j) / chr.sequenceLen * 2 * pi;
                largeArc = dTheta >= pi;
                r = r0(iseven(posStrnds(j, 2)) + 1);
                if posStrnds(j, 2) <= 2
                    cx = centerX(1);
                    cy = centerY(1);
                else
                    cx = centerX(2);
                    cy = centerY(2);
                end
                
                strokeColor = strokeColor_set;
                if ~isequal(helicase,[0;0]) && helicase(2) > 0
                    if posStrnds(j,2) == 2 && (posStrnds(j,1) > helicase(1) || posStrnds(j,1) < helicase(2))
                        strokeColor = 'red';       
                    elseif posStrnds(j,2) == 3 && (posStrnds(j,1) > helicase(1) || posStrnds(j,1) < helicase(2))
                          strokeColor = 'red';       
                    end
                end
                if lens(j) == chr.sequenceLen
                    svg = [svg sprintf('<circle cx="%d" cy="%d" r="%d" stroke="%s" stroke-width="%f" fill="none"/>\n', ...
                        cx, cy, r, ...
                        strokeColor, strokeWidth)]; %#ok<*AGROW>
                else
                    svg = [svg sprintf('<path d="M%f,%f A%f,%f %d %d,%d %f %f" stroke="%s" stroke-width="%d" fill="none"/>\n', ...
                        cx + r*cos(theta), cy + r*sin(theta), ...
                        r, r, 0, largeArc, 1, cx + r*cos(theta+dTheta), cy + r*sin(theta+dTheta), ...
                        strokeColor, strokeWidth)];
                end
            end
        end
        
        function svg = drawCircles(this, posStrnds, lens, r0, strokeColor, strokeWidth, helicase, centerX, centerY, rad, fill)
           chr = this.simulation.state('Chromosome');
%             strokeColor_set = strokeColor;
            
            svg = [];
            for j = 1:size(posStrnds, 1)
                theta = posStrnds(j, 1) / chr.sequenceLen * 2 * pi - pi/2;
                dTheta = lens(j) / chr.sequenceLen * 2 * pi;
                largeArc = dTheta >= pi;
                r = r0(iseven(posStrnds(j, 2)) + 1);
                if posStrnds(j, 2) <= 2
                    cx = centerX(1) + r*cos(theta);
                    cy = centerY(1) + r*sin(theta);
                else
                    cx = centerX(2) + r*cos(theta);
                    cy = centerY(2) + r*sin(theta);
                end
                
%                 strokeColor = strokeColor_set;
%                 if ~isequal(helicase,[0;0]) && helicase(2) > 0
%                     if posStrnds(j,2) == 2 && (posStrnds(j,1) > helicase(1) || posStrnds(j,1) < helicase(2))
%                         strokeColor = 'red';       
%                     elseif posStrnds(j,2) == 3 && (posStrnds(j,1) > helicase(1) || posStrnds(j,1) < helicase(2))
%                           strokeColor = 'red';       
%                     end
%                 end

                svg = [svg sprintf('<circle cx="%d" cy="%d" r="%d" stroke="%s" stroke-width="%f" fill="%s"/>\n', ...
                    cx, cy, rad, ...
                    strokeColor, strokeWidth, fill)]; %#ok<*AGROW>
            end
        end
        
        function renderFrames(this)
            switch this.frameRenderer
                case 'librsvg'
                    this.renderFrames_librsvg();
                case 'imagemagick'
                    this.renderFrames_imagemagick();
                case 'inkscape'
                    this.renderFrames_inkscape();
                case 'batik'
                    this.renderFrames_batik();
                otherwise
                    throw(MException('CellGeometryAnimation:error', 'Undefined renderer %s', this.frameRenderer))
            end
        end
        
        function renderFrames_librsvg(this)
            tic
            nFrames = floor(numel(this.states.Time.values) / this.downsampleStepSizeSec);
            renderFileName = [this.tmpDirectory, filesep 'renderFrames.sh'];
            fid = fopen(renderFileName, 'w');
            fprintf(fid, '#!/bin/sh\n');
            for frameIdx = 1:nFrames
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
            toc
        end
        
        function renderFrames_imagemagick(this)
            [status, result] = system(sprintf([...
                'mogrify '...
                '-format png -define png:bit-depth=8 -define png:color-type=2 -define png:compression-level=0 -define png:compression-strategy=0 '...
                '-density %d -resize %dx%d '...
                '"%s%s%s.svg"'], ...
                round(3*this.movieDPI/90*72), round(this.movieWidth), round(this.movieHeight), this.tmpDirectory, filesep, strrep(this.frameFilePattern, '%d', '*')));
            if status ~= 0
                throw(MException('Animation:error', 'Failed to convert svg to png: %s', result));
            end
        end
        
        function renderFrames_inkscape(this)
            if(isunix())
                a=getenv('LD_LIBRARY_PATH');
                a=['/usr/lib:' a];
                setenv('LD_LIBRARY_PATH',a);
            end
            nFrames = floor(numel(this.states.Time.values) / this.downsampleStepSizeSec);
            for frameIdx = 1:nFrames
                %print status
                if this.verbosity > 0
                    fprintf('Drawing frame %d of %d frames\n', ...
                        frameIdx, nFrames);
                end
                [status, result] = system(sprintf('inkscape "%s" --export-png="%s" --export-area-page -w %d -h %d', ...
                    sprintf(['%s%s' this.frameFilePattern '.svg'], this.tmpDirectory, filesep, frameIdx), ...
                    sprintf(['%s%s' this.frameFilePattern '.png'], this.tmpDirectory, filesep, frameIdx), ...
                    round(this.movieWidth), round(this.movieHeight)));
                if status ~= 0
                    throw(MException('Animation:error', 'Failed to convert svg to png: %s', result));
                end
            end
        end
        
        function renderFrames_batik(this)
            [status, result] = system(sprintf('java -jar batik-rasterizer.jar -m image/png -d "%s" "%s%s%s.svg"', ...
                this.tmpDirectory, this.tmpDirectory, filesep, strrep(this.frameFilePattern, '%d', '*')));
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
            fprintf(this.subtitleFid, '%s\n\n', datestr(this.states.Time.values(timeIdx) / 3600 / 24, 'HH:MM:SS'));
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
                    throw(MException('CellGeometryAnimation:error', 'Undefined renderer %s', this.movieRenderer))
            end
        end
        
        function renderMovie_ffmpeg(this)
            cmd = sprintf([...
                'ffmpeg -f image2 -i "%s/%s.png" ' ...
                '-vcodec %s -r %f -b %dk -loop_output %d ' ...
                '-an ' ...
                '-timestamp now -metadata title="%s" -metadata description="%s" -metadata author="%s" -metadata year="%s"  -metadata copyright="%s"' ...
                '-y "%s"'
                ], ...
                this.tmpDirectory, this.frameFilePattern, ...
                this.codec.ffmpeg, this.frameRate, this.bitRate, this.loopCount, ...
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
                'mencoder "mf://%s.png" -mf type=png:fps=%d -sub "%s" -o "%s" ' ...
                '-ovc lavc -lavcopts vcodec=%s:vbitrate=%d:vpass=1 ' ...
                '-subcp enca:en:iso-8859-7 ' ...
                '-info name="%s":subject="%s":artist="%s":copyright="%s" ' ...
                ], ...
                this.tmpDirectory, this.frameFilePattern, this.frameRate, this.subtitleFileName, movieFileName, ...
                this.codec.mencoder, this.bitRate, ...
                this.title, this.description, this.author, this.copyright);
            [status, result] = system(cmd);
            if status ~= 0
                throw(MException('Simulation:error', 'Failed to render movie: %s', result));
            end
        end
    end
end