function pz_plot_dcm(Ep,xy,names,width_range)
% Plots a DCM connectivity matrix
%
% Ep          - [nxn] connectivity matrix for n regions
% xy          - [nx2] coordinate matrix of locations for each region. These 
%               will be rescaled to the range [-1.5 1.5]
% names       - {1xn} name of each region (optional)
% width_range - [1 x 2] matrix with the lower and upper range of widths to
%               use (optional, default [5 20])
%
% Example:
% figure('Color','w');
% xy    = [DCM.xY.xyz]';
% names = {DCM.xY.name};
% pz_plot_dcm(DCM.Ep.A,xy,names);

% ----
n = size(Ep,1);

% Prepare names
if nargin < 3 || isempty(names)
    names = cell(1,n);
    for i = 1:n
        names{i} = sprintf('R%d',i);
    end
end

assert(size(xy,1)==n,'xy needs to have one row per region');
assert(length(names)==n,'names should contain one element per region');

% Remove 3rd coordinate dimension if provided
if size(xy,2) > 2
    xy = xy(:,1:2);
end

% Rescale coordinates
xy(:,1) = pz_rescale(xy(:,1), -1.5, 1.5);
xy(:,2) = pz_rescale(xy(:,2), -1.5, 1.5);

% Remove v. small parameters
is_zero = (abs(Ep) < 0.001);

% Compute widths
if nargin < 4
    width_range = [5 20];
end
widths = Ep;
widths(~is_zero) = ...
    pz_rescale(abs(widths(~is_zero)),width_range(1),width_range(2));

if length(unique(Ep(~is_zero))) < 2
    widths(~is_zero) = 2;
end

% Prepare axes
% -------------------------------------------------------------------------
axis(gca);
axis equal;
axis off;

% Plot connections
% -------------------------------------------------------------------------
for r_to = 1:n
    for r_from = 1:n
        
        A = xy(r_from,:);
        B = xy(r_to,:);
        
        if abs(Ep(r_to,r_from)) > 0.001
            if Ep(r_to,r_from) > 0
                colour = [45, 164, 33]/255;
            else
                colour = [200, 56, 41]/255;
            end
            arrowsize = widths(r_to,r_from);
            
            if r_to == r_from 
                curve = -0.5;
                B = B + [0.0001 0.0001];
            else
                curve = 0.1;
            end
            draw_curved_arrow(A,B,curve,colour,arrowsize);            
        else
            % nothing
        end
    end
end

% Plot labels
% -------------------------------------------------------------------------
for i = 1:n
    h  = circle(xy(i,1), xy(i,2),0.5,'FaceColor','k'); hold on;
    h2 = text(xy(i,1), xy(i,2), names{i},'Color','w',...
        'HorizontalAlignment','center','FontWeight','Bold');    
end

% Update axis limits
% -------------------------------------------------------------------------
hold off;

limx = [min(xy(:,1)) max(xy(:,1))];
limx = limx + [-0.5 0.5];

limy = [min(xy(:,2)) max(xy(:,2))];
limy = limy + [-0.5 0.5];

xlim(limx);
ylim(limy);

% -------------------------------------------------------------------------
function h= circle(x,y,r,varargin)
h = rectangle('Position',[x-(r/2) y-(r/2) r r],'Curvature',[1 1],varargin{:});

% -------------------------------------------------------------------------
function draw_curved_arrow(A,B,curviness,color,size,varargin)
% Thanks to: 
% http://stackoverflow.com/questions/27460081/how-do-i-know-position-of-point-perpendicular-to-line

% Vector from A->B
AB = B - A;

% Midpoint
M = A + AB / 2;

% Find a vector orthogonal to AB and to an arbitrary vector v
v = [0 0 1];
V = cross([AB 0], v);

% Set length of orthogonal vector
V = (V / norm(V)) * curviness;

% Travel along orthogonal vector
c = [M 0] + V;

% Draw curve
P = bezier(A, c(1:2), B, 'Color', color, varargin{:});

% Find midpoint of bezier
np = length(P);
MB = P(ceil(np * 0.8),:);
MB2 = P(ceil(np * 0.8+1),:);

% Put arrow at the midpoint pointing to destination
h=arrowh([MB(1) MB2(1)],[MB(2) MB2(2)],color,size*100);

% -------------------------------------------------------------------------
function new_X = pz_rescale(X, new_min, new_max)

current_min = min(X(:));
current_max = max(X(:));

scale_factor = (current_max - current_min) / (new_max - new_min);
new_X = new_min + (X - current_min) / scale_factor;    