function [circData, circHandles, frame, S] = bubblebath(S)
% Create a 2D plot of circles (or 'bubbles') with random centers and varying radii. 
% BUBBLEBATH() uses default parameter values.
% BUBBLEBATH(S) uses a structure S with the following fields, all of which are optional. 
%     axisHandle: axis handle; an axis is created in a new figure if missing.
%      frameSize: 1x2 vector defining the frame's [width,height] (centered at 0).   default [50,50]
%       circSize: 1x2 vector defining [min, max] range of radii.                    default [0.2,5]
%                 1xn vector defining all possible radii when nSizes is NaN.
%         nSizes: Numeric scalar, number of radii to include within circSize.       default 25
%                 NaN: indicates that circSize defines all radii.
% maxCircsPerRad: Maximum number of circles per radius to avoid memory problems.    default 5000
%                 Positive integer or inf.
%     circPoints: Scalar integer greater or equal to 3 describing the number of     default 628
%                 coordinates used to generate each polygon between [0, 2*pi].
%                 Use S.circPoints=3 to draw triangles, = 4 for squares, etc.
%          maxIt: Numeric scalar, max number of attempts to draw circle per radius. default 200
%       edgeType: 0:circle edges can expand outside the frame.                      default 2
%                 1:circles must be entirely inside the frame.
%                 2:circle edges that would expand beyond the frame are cut off.
%                 3:circles that intersect the frame are wrapped to the other side.
%        density: Density of circles; numeric scalar [0 < density <= 1].            default 0.7
%        overlap: True: circles can overlap any amount.                             default false
%                 False|0: circle edges can touch but not overlap.
%                 Positive scalar: sets the minimum edge distance (see overlapType).
%                 Negative scalar: sets the maximum edge overlap (see overlapType).
%    overlapType: 'relative': overlap amount is r*overlap, r is radius of cirlce.   default 'relative'
%                 'absolute': overlap amount is constant for all radii.
% supressWarning: true: supress internal warnings; false: show warning.             default false
%      drawFrame: true: draw a rectangle showing the frame; false: don't draw it.   default true
%
% Outputs
%   * circData: [m x 3] matrix of m circles data showing [xCenter, yCenter, radius].
%   * circHandles: [m x 1] vector of handles for each line object / circle.
%   * frame: handle to 'frame' rectangle (GraphicsPlaceHolder if drawFrame is false).
%   * S: A structure with fields listing all parameters used to reproduce the figure
%       also including S.rng which is the random number generator state you can use
%       to reproduce a figure.
%
% Examples
%   For many examples, see <a href = "https://www.mathworks.com/matlabcentral/fileexchange/70348">bubblebath_examples.mlx</a>
%   Example 1: draw bubblebath plot using default inputs.
%       BUBBLEBATH()
%   Example 2: specify some of the inputs
%       S.circSize = [0.5, logspace(0, 1, 5)];
%       S.nSizes = NaN;
%       BUBBLEBATH(S)
%   Example 3: Reproduce an exact replica
%       [~,~,~,Sout] = bubblebath();
%       newFig = figure();
%       Sout.axisHandle = axes(newFig);
%       rng(Sout.rng) % seed the random number generator
%       bubblebath(Sout)
%
% Examples of how to use outputs
%   * Calculate the area for each circle: circArea = pi * circData(:,3).^2;
%   * Change color and width of circle edges:  set(circHandles, 'color', 'r', 'LineWidth', 3)
%   * Fill circles with color:
%         c = jet(size(circData,1)); %choose colors
%         for i = 1:size(circData,1)
%             rectangle('Position', [circData(i,1:2)-circData(i,3), circData(i,[3,3])*2], 'Curvature',[1,1], 'FaceColor',c(i,:))
%         end
%   * Count number of circles for each circle size:
%       [g, radius] = findgroups(circData(:,3));
%       count = splitapply(@length, circData(:,3), g);
%       table(radius, count)
%
% Details
%   * For each radius, the number of circles we attempt to draw is a proportion between
%       the circle's area and the frame's area. For 'n' circle sizes, the total area of
%       the frame is evenly split into 'n' portions and the number of circles per circle
%       size is merely framePortion / circleArea, rounded.
%   * Circles are drawn in order of circle size in descending order. The distance
%       between circles is calculated by the distance between all circle centers minus
%       the radii.  That gives the distance between all circle edges.
%   * Since circle centers are chosen randomly, we check that the edges won't overlap
%       with existing circles.  If they do, we choose another set of random centers
%       and that process repeats in a while-loop until either all circle centers
%       do not overlap or we meet a maximum iteration (set by user).  A warning
%       message appears if the while loop expired without drawing all expected
%       circles.
%
% Requires Matlab r2016a or later.
% Copyright (c) 2019, Adam Danz
% All rights reserved
% Source: <a href = "https://www.mathworks.com/matlabcentral/fileexchange/70348">bubblebath</a>
% Author: <a href = "https://www.mathworks.com/matlabcentral/profile/authors/3753776-adam-danz">Adam Danz</a>
% This function was written for a Matlab forum question:
% https://www.mathworks.com/matlabcentral/answers/446114-non-overlapping-random-circles
% Revision history
% vs 1.0.0  02/21/2019  Initial FEX upload
% vs 1.2.0  09/02/2019  Added waitbar, adapted to work with r2016b.
% vs 2.0.0  02/02/2020  Converted to proper function with input options; added several new
%                       options: supressWarning, overlap, overlapType, maxCircsPerRadius,
%                       axisHandle; names of some variables; input validation; documentation.
% vs 2.1.0  02/03/2020  S.nSizes=nan was being overwritten which affected the S output. Fixed.
%                       Updated waitbar message in include max number of circles being drawn.
% vs 2.2.0  04/07/2020  S.circPoints added (ability to draw n-sided polygons). Unrecognized
%                       fields of S will now throw warning rather than an error. fs2 renamed
%                       intlFrame and changed to 2xn mat so S.edgeType=1 fits more bubbles.
%                       Max iterations warning now without backtrace. maxCircsPerRadiuschanged
%                       to maxCircsPerRad but both are accepted for backward compat.
% vs 2.3.0  01/21/2021  Added edgeType 3; added approveNewCircle() to clean up while-loop.
%                       Algo change: replacing only the random circs that aren't already isOK.
%                       This change affects the reproducibility of plots from older versions
%                       with the same rng seed. 'd' renamed to 'nCirc' and redefined to consider
%                       intlFrame and minDist [1]. Added bubblebath_examples.mlx.
%% Default inputs
% Fill in missing params with default values.
if nargin == 0
    S = struct();
end
if ~isfield(S,'axisHandle') || isempty(S.axisHandle)
    % Create figure & axes
    fh = figure();
    S.axisHandle = axes(fh); % r2016a or later
end
if ~isfield(S,'frameSize') || isempty(S.frameSize)
    % Frame size (arbitrary units): size of axes, centered at (0,0).
    S.frameSize = [50, 50]; %[width, height]
end
if ~isfield(S,'circSize') || isempty(S.circSize)
    % Circle sizes: Select the minimum and maximum radii
    % Radii are linearly spaced. ie: linspace(S.circSize(1), S.circSize(2), S.nSizes).
    % Alternatively, specify all radii directly in a 1xn vector (nSizes must be NaN).
    S.circSize = [.2, 5]; %[smallest, largest] radius
end
if ~isfield(S,'nSizes') || isempty(S.nSizes)
    %number of circle sizes between S.circSize(1) and S.circSize(2)
    % When NaN, circSize specifies radii directly.
    S.nSizes = 25;
end
if isfield(S, 'maxCircsPerRadius')
    % vs 2.0 mistakenly used maxCircsPerRad in help section but the
    % true field name was maxCircsPerRadius.  From vs 2.2 we'll
    % accept both and will change fieldname here.
    S.maxCircsPerRad = S.maxCircsPerRadius;
    S = rmfield(S, 'maxCircsPerRadius');
end
if ~isfield(S,'maxCircsPerRad') || isempty(S.maxCircsPerRad)
    % Set the maximum number of circles per radius.  The squareform
    % function may cause memory problems as this number gets larger
    % (eg ~12000 for 64bit Windows 10).  Use inf for unlimited.
    S.maxCircsPerRad = 5000;
end
if ~isfield(S,'circPoints') || isempty(S.circPoints)
    % Set the number of coordinates used to draw each circle [0,2*pi].
    S.circPoints = 628;
end
if ~isfield(S,'maxIt') || isempty(S.maxIt)
    % max iterations: how many attempts should be made to find circle
    % locations that do not overlap with other circles?
    S.maxIt = 200;
end
if ~isfield(S,'edgeType') || isempty(S.edgeType)
    % Decide whether circles shoule be entirely inside of the frame
    % 0 = Cirlces edges can expand outside of the frame
    % 1 = Cirlces must be entirely inside of the frame
    % 2 = Circle edges that extend beyond the frame are cut off
    % 3 = Circles at edges are wrapped
    S.edgeType = 2;
end
if ~isfield(S,'density') || isempty(S.density)
    % Density of circles:
    % a value greater than 0 but less than or equal to 1.
    % As the value approaches 0, less circles will be drawn.
    S.density = .7;
end
if ~isfield(S,'overlap') || isempty(S.overlap)
    % Allow overlap: True/False or a numeric scalar.
    % When true, bubbles can overlap any ammount.
    % When false or 0, bubble edges can touch but cannot overlap.
    % When numeric & S.overlapType is 'absolute', then all bubble edges will
    % have a minimum edge distance of S.overlap. When S.overlabType
    % is 'relative', then all bubble edges will have a min edge distance
    % of r*S.overlap where r is the radius of the smaller circle
    % being added. Negative values specify the amount of overlap allowed.
    S.overlap = false;
end
if ~isfield(S,'overlapType') || isempty(S.overlapType)
    % See description for overlap.  This parameter is
    % ignored when overlap is 0, true, or false.
    S.overlapType = 'relative';
end
if ~isfield(S,'supressWarning') || isempty(S.supressWarning)
    % When true, the internal warning will be supressed.
    S.supressWarning = false;
end
if ~isfield(S,'drawFrame') || isempty(S.drawFrame)
    % Draw the defined frame around the bubble plot.
    S.drawFrame = true;
end
% input validation
assert(ishghandle(S.axisHandle,'axes'),'axisHandle must be an axis handle.')
validateattributes(S.frameSize,{'double'},{'size',[1,2],'>',0},mfilename,'frameSize')
assert(isnumeric(S.nSizes) && isscalar(S.nSizes) && (isnan(S.nSizes) || mod(S.nSizes,1)==0),...
    'nSizes is expected to be a scalar, nonzero, positive integer or NaN.')
if isnan(S.nSizes)
    assert(isrow(S.circSize) && isnumeric(S.circSize) && all(~isinf(S.circSize)) && all(S.circSize>0), ...
        'When nSizes is NaN, circSize is expected to be a 1xn numeric vector, with finite, positive, nonzero values.')
elseif isnumeric(S.circSize) && numel(S.circSize) > 2
    error('When defining all possible radii in S.circSize, S.nSizes must be NaN.')
else
    validateattributes(S.circSize,{'double'},{'size',[1,2],'>',0,'increasing'},mfilename,'circSize')
end
assert(isnumeric(S.maxCircsPerRad) && isscalar(S.maxCircsPerRad) && S.maxCircsPerRad>0 ...
    && (mod(S.maxCircsPerRad,1)==0 || isinf(S.maxCircsPerRad)), 'maxCircsPerRad must be a positive, non-zero integer or inf.')
validateattributes(S.circPoints,{'double'},{'numel',1,'>=',3,'integer'},mfilename,'circPoints')
validateattributes(S.maxIt,{'double'},{'numel',1,'>',0,'integer'},mfilename,'maxIt')
assert(ismember(S.edgeType,0:3),'edgeType must be 0, 1, 2, or 3.')
validateattributes(S.density,{'double'},{'numel',1,'>',0,'<=',1},mfilename,'density')
assert((islogical(S.overlap)||isnumeric(S.overlap))&&numel(S.overlap),'overlap must be true, false, or a numeric scalar value.')
assert(any(strcmpi(S.overlapType,{'absolute','relative'})),'overlapType must be either ''relative'' or ''absolute''.')
validateattributes(S.supressWarning,{'logical'},{'numel',1},mfilename,'supressWarning')
validateattributes(S.drawFrame,{'logical'},{'numel',1},mfilename,'drawFrame')
% Check for included parameters (fields of S) that are not recognized.
acceptedFields = {'axisHandle','frameSize','circSize','nSizes','circPoints','maxIt','edgeType','density',...
    'overlap','overlapType','supressWarning','drawFrame','maxCircsPerRad'};
params = fields(S);
foreignFields = ismember(params,acceptedFields);
if ~all(foreignFields) && ~S.supressWarning
    warning('Parameter(s) [%s] not recognized and will be ignored.', strjoin(params(~foreignFields),', '))
end
% Store rng state for reproducibility
S.rng = rng();
%% Add soap
hold(S.axisHandle,'on')
axis(S.axisHandle,'equal') %set aspect ratio to 1:1
xlim(S.axisHandle, S.frameSize(1)/2 * [-1.05, 1.05]) %with some extra space
ylim(S.axisHandle, S.frameSize(2)/2 * [-1.05, 1.05]) %with some extra space
% determine minimum distance between circles allowed
if islogical(S.overlap) && ~S.overlap
    minDist = 0;
elseif islogical(S.overlap) && S.overlap
    minDist = -inf;
else
    minDist = S.overlap;
end
% determine circle sizes (largest to smallest)
if isnan(S.nSizes)
    nSizes = numel(S.circSize);
    r = sort(S.circSize,'descend');
else
    nSizes = S.nSizes;
    r = linspace(S.circSize(2), S.circSize(1), nSizes);
end
% Identify the internal frame of possible circle centers (intlFrame)
switch S.edgeType
    case 0  % Cirlce edges can expand outside of the frame
        intlFrame = repmat(S.frameSize(:), 1, nSizes);
        clearBorder = 0;
    case 1  % Circle must be entirely inside of the frame
        intlFrame = bsxfun(@minus, S.frameSize(:), r*2);
        clearBorder = 0;
    case 2  % Circle edges end at the frame
        intlFrame = repmat(S.frameSize(:), 1, nSizes);
        clearBorder = 1;
    case 3 % Circle edges are wrapped
        intlFrame = repmat(S.frameSize(:), 1, nSizes);
        clearBorder = 2;
end
% adjust min distance between circles based on overlap type
if strcmpi(S.overlapType,'absolute')
    minDist = repmat(minDist,size(r));
elseif strcmpi(S.overlapType,'relative')
    minDist = r*minDist;
end
% determine approximate number of circles to draw for each size (see [1])
nCirc = ceil(prod(intlFrame,1)./nSizes ./ (pi * r.^2) * S.density);
if any(nCirc > S.maxCircsPerRad) && ~S.supressWarning
    warning('The maximum number of circles for radii [%s] have been limited to %d.', ...
        regexprep(strtrim(num2str(r(nCirc>S.maxCircsPerRad))),' +',' '), S.maxCircsPerRad)
end
nCirc(nCirc > S.maxCircsPerRad) = S.maxCircsPerRad;
% Throw error if largest circle is larger than frame
frameArea = prod(S.frameSize);
circAreas = pi * r.^2;
assert(max(circAreas) <= frameArea, 'The area of the largest circle (%.1f) is larger than the frame''s area (%.1f)', max(circAreas), frameArea)
% Loop through each circle size
circdata = []; %[xCenter, yCenter, radius] of each drawn circle
h = cell(nSizes,1); % handles to the line objects for each circle
wb = waitbar(0,'initializing...','name',mfilename);
originalWarnState = warning('backtrace');
warning backtrace off
for i = 1:nSizes
    % Reset circle & iteration count
    cCount = 0;
    iCount = 0;
    xRand = zeros(nCirc(i),1);
    yRand = xRand;
    isOK = false(size(xRand));
    % Keep drawing random coordinates until either all circs
    % are drawn or we reach max number of attempts
    while cCount < nCirc(i) && iCount < S.maxIt
        % Randomly choose center coordinates
        xRand(~isOK) = (rand(nCirc(i)-cCount,1) - 0.5) * intlFrame(1,i);
        yRand(~isOK) = (rand(nCirc(i)-cCount,1) - 0.5) * intlFrame(2,i);
        xyr = [xRand, yRand, repmat(r(i), size(yRand))];
        
        % For edge-wrapping
        if clearBorder==2
            [xRand, yRand, isOK, cCount] = wrapEdges(circdata, xyr, minDist(i), nCirc(i), isOK, S);
        else
            % Determine if new circles overlap with others or themselves
            [isOK, cCount, xRand, yRand] = approveNewCircle(circdata, xyr, minDist(i), true);
        end
        
        % Peak at random circles (for development purposes only)
        %   peakAtCircle([xRand,yRand,repmat(r(i), size(yRand))],S)
        
        iCount = iCount + 1; %iteration count
        % Update waitbar
        if ishghandle(wb)
            waitbar(max(iCount/S.maxIt,cCount >= nCirc(i)),wb,sprintf('Trying to find space for up to %d circles with radius = %.2f\nFinal radius: %.2f',nCirc(i),r(i),r(end)));
        end
    end
    % If we had to quit searching, throw warning.
    if iCount >= S.maxIt && ~S.supressWarning
        warning('Max iteration reached. %d of the requested %d circles drawn for radius %.3f', cCount, nCirc(i), r(i))
    end
    % Store all final circle data
    circdata = [circdata; [xRand(isOK), yRand(isOK), repmat(r(i), sum(isOK), 1)]]; %#ok<AGROW>
    % Draw circles
    if any(isOK)
        h{i} = drawcircles([xRand(isOK), yRand(isOK), repmat(r(i), sum(isOK), 1)], clearBorder, S);
    end
    
end
warning(originalWarnState) %return original warning state
% Draw frame
if S.drawFrame
    frame = rectangle(S.axisHandle, 'position', [-S.frameSize/2, S.frameSize], 'LineWidth', 2);
else
    frame = gobjects(1);
end
circHandles = [h{:}]';
% Remove waitbar
if ishghandle(wb)
    delete(wb)
end
if nargout > 0
    % Produce output only when requested
    circData = circdata;
end
%% Helper functions
function peakAtCircle(xyr,S) %#ok<DEFNU>
% For development and troubleshooting only; creates a rough plot of circles and labels them.
% xyr is mx3, [x,y,radius]; S is input struct.
figure('Name',mfilename)
hold on
nanIdx = any(isnan(xyr),2);
rectangle('Position',[-S.frameSize(1)/2,-S.frameSize(2)/2,S.frameSize],'LineWidth',2)
arrayfun(@(x,y,r)rectangle('Position',...
    [x-r,y-r,r.*[2,2]],'Curvature',[1,1],'EdgeColor','r',...
    'LineStyle','--'),xyr(~nanIdx,1),xyr(~nanIdx,2),xyr(~nanIdx,3));
text(xyr(:,1),xyr(:,2),compose('%d',1:size(xyr,1)),'Vert','Middle','Horiz','Center','FontSize',8)
axis equal; grid on;
xlim(S.frameSize./[-1.6,1.6]);
ylim(xlim);
set(gca,'xtick',get(gca,'ytick'))
title(sprintf('%s.m test & development', mfilename))
function [x, y, isOK, cCount] = wrapEdges(circdata, xyr, minDist, minCount, doNotCopy, S)
% Detects circles that cross the frame and produces a second set of circles
% that are wrapped to the other x/y axis.  Then tests for approval and returns
% the accepted xyr values.
%  circdata: nx3 defs for existing circles [x,y,radius]
%  xyr: mx3 defs for new circles to be tested [x,y,radius]
%  minDist: scalar, min distance appepted between circles
%  minCount: minimum number of circles to return, even if they are NaNs.
%  doNotCopy = logical col vec with height size(xyr,1) ID'ing which circs are already OK and don't copy.
%  S: the input structure to the main fcn.
% Peak at random circles (for development purposes only)
%   peakAtCircle(xyr,S)
% isXYCrossing is % nx2 for n circles, [x,y]; needs to support r2016a
isXYCrossing = abs(xyr(:,[1,2])) > bsxfun(@minus, S.frameSize/2, xyr(:,3)) & ~doNotCopy*[1,1]; 
% Create a second or third set of circles for those that are border crossers
if any(isXYCrossing(:))
    % Compute distance-to-edge and new wrapped centers for *all* circles (even if they don't need wrapped)
    distanceToBorder = bsxfun(@minus,S.frameSize/2, abs(xyr(:,1:2))).*sign(xyr(:,1:2)); % signed distance, neg dist to right/bottom border
    newCenters = bsxfun(@times,-S.frameSize/2,sign(xyr(:,1:2))) - distanceToBorder;
    % Create centers for new circles; x and y must be replicated independently.
    % If a circle breaks 1 edge, another circle will be defined.  If a circle
    % breaks 2 edges, two circles will be defined.
    xyrAdded = [
        [newCenters(isXYCrossing(:,1),1), xyr(isXYCrossing(:,1),2:end)];
        [xyr(isXYCrossing(:,2),1), newCenters(isXYCrossing(:,2),2), xyr(isXYCrossing(:,2),3)]
        ];
    [rowCopied, ~] = find(isXYCrossing); % xyrAdded(j,:) is a copy of xyr(rowCopied(j),:)
    % update list of circles to test
    isCopied = [false(size(xyr,1),1); true(size(xyrAdded,1),1)];
    xyr = [xyr; xyrAdded];
else
    isCopied = false(size(xyr,1),1);
    rowCopied = nan(size(isXYCrossing,1),1);
end
% Test for approval
isOK = approveNewCircle(circdata, xyr, minDist, false); %sorted later
% If any wrapped circles are not accepted, eliminate their original counterpart
isOK(rowCopied(~isOK(isCopied))) = false;
% If any original circles were not accepted, eliminate their wrapped counterparts
wrappedRows = find(isCopied);
isOK(wrappedRows(ismember(rowCopied, find(~isOK(~isCopied))))) = false;
% Circle count only inlude the original portion, not the wrapped portion.
cCount = sum(isOK(~isCopied));
% flag unaccepted circles for removal using NaN, sort so accepted on top,
% make sure there's enough trailing space holders for the number of remaining
% circles needed.
nRemaining = max(minCount-cCount,0);
x = [xyr(isOK,1); NaN(nRemaining,1)];
y = [xyr(isOK,2); NaN(nRemaining,1)];
isOK = [true(sum(isOK),1); false(nRemaining,1)];
function [isOK, cCount, xRand, yRand] = approveNewCircle(circdata, xyr, minDist, returnSorted)
% Determine if new circles overlap with others or themselves.
% INPUTS
%   circdata: nx3 defs for existing circles [x,y,radius]
%   xyr: mx3 defs for new circles to be tested [x,y,radius]
%   minDist: scalar, min distance appepted between circles
%   returnSorted: scalar logical; only wrapEdges()>approveNewCircle() should use false!
% OUTPUTS:
%   isOK: mx1 logical vector accepting/denying the *sorted* xyr circles
%       where accepted circles are moved to the top when returnSorted is true
%   cCount: scalar integer, number of accepted circles
%   xRand,yRand: the circ centers in same order as isOK (returnSorted must be true).
% Peak at random circles (for development purposes only)
%   peakAtCircle(xyr,S)
allCircs = [circdata; xyr];
radMat = repmat(allCircs(:,3),1,size(allCircs,1)) + repmat(allCircs(:,3).',size(allCircs,1),1); %changed 190901 to work with r2016a
dist = tril(squareform(pdist(allCircs(:,1:2))) - radMat, -1);
if isempty(dist)
    isOK = true; %when xyr only has 1 row
else
    dist(triu(true(size(dist)))) = inf;
    isOK = all(dist(size(circdata,1)+1:end,:) >= minDist, 2);
    %Put accepted coordinates on top so only unaccpted circles are replaced
    if returnSorted
        xyr = [xyr(isOK,:); xyr(~isOK,:)];
        isOK = sort(isOK,'descend');
    end
end
cCount = sum(isOK);  %cirlce count for current radius
xRand = xyr(:,1);
yRand = xyr(:,2);
function h = drawcircles(xyr, clearBorder, S)
% Draw circle given center and radius
ang = linspace(0, 2*pi, S.circPoints+1);
xp = xyr(:,3)*cos(ang) + repmat(xyr(:,1),1,numel(ang)); %changed 190901 to work with r2016a
yp = xyr(:,3)*sin(ang) + repmat(xyr(:,2),1,numel(ang)); %changed 190901 to work with r2016a
if clearBorder==1 || clearBorder==2
    % remove data outside of frame, if requested
    xp(abs(xp) > S.frameSize(1)/2) = NaN;
    yp(abs(yp) > S.frameSize(2)/2) = NaN;
end
h = plot(S.axisHandle, xp', yp', 'k')';
%% Footnotes
% [1] nCirc used to be named 'd' and was defined by (2:nSizes+1).^2; Then prior
%   to vs. 2.3.0 it was defined by the line below and renamed to nCirc.
%   nCirc = ceil(ceil((frameArea / nSizes) ./ circAreas) * S.density); 