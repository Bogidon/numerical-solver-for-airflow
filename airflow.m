%% Linarity II Final: Numerical Solver for Airflow
% Authors: Bogdan Vitoc, Christopher Lather, William Fairman
% Date: December 11th, 2018
% Homepage: https://github.com/Bogidon/numerical-solver-for-airflow
%
% airflow.m: a program that calculates 2D airflow around any
% inputs:
%   - spacing: distance between elements on the computation grid
%   - left_boundary_potential: velocity potential at the left boundary that
%   is held constant.
%   - grid and shape parameters:
%
%           shape_xmin  shape_xmax               
%                ?           ?                   
%              ?????????????????????  grid_ymax  
%              ? ?    grid   ?   ?               
%              ? ?           ?   ?               
% shape_ymax ?????????????????   ?               
%              ? ?shape      ?   ?               
%              ? ?(shapepath)?   ?               
% shape_ymin ?????????????????   ?               
%              ?                 ?               
%              ????????????????????? grid_ymin   
%              ?                 ?               
%          grid_xmin         grid_xmax                                                          

%% figures
fig1 = figure;
fig2 = figure;
fig1.WindowStyle = 'docked';
fig2.WindowStyle = 'docked';

%% simulation (required for the other sections)
% input variables
spacing = 0.1;
grid_xmin = -4;
grid_xmax = 4;
grid_ymin = -20;
grid_ymax = 20;
shape_xmin = -1.5;
shape_xmax = 4.5;
shape_ymin = -2.25;
shape_ymax = 2.25;
imagepath = 'shapes/swallow.jpg';
left_boundary_potential = 50;

% create the scalar grid
x = grid_xmin:spacing:grid_xmax;
y = grid_ymin:spacing:grid_ymax;
[xx,yy] = meshgrid(x,y);
[exterior, ghosts, ghosts_struct, interior] = elliptic_curve(shape_xmin,shape_xmax,shape_ymin,shape_ymax,imagepath,x,y);
left_side = logical(padarray(ones(length(y),1),[0,length(x)-1],'post'));
right_side = logical(padarray(ones(length(y),1),[0,length(x)-1],'pre'));
top_side = logical(padarray(ones(1,length(x)),[length(y)-1,0],'post'));
bottom_side = logical(padarray(ones(1,length(x)),[length(y)-1,0],'pre'));
sides = left_side | right_side |top_side | bottom_side;
exterior = xor(exterior,sides);
Z = zeros(size(exterior)) + (left_side * left_boundary_potential);

for i = 1:10000
    Z_old = Z;
    
    % relaxation
    for jx = 1:size(Z,1)
        for jy = 1:size(Z,2)
            if exterior(jx,jy)
                Z(jx,jy) = (Z_old(jx-1,jy) + Z_old(jx+1,jy) + Z_old(jx,jy-1) + Z_old(jx,jy+1)) / 4;
            end
        end
    end
    
    % ghost points = mirror
    for ghost = ghosts_struct
        Z(ghost.row,ghost.col) = mean(Z_old(ghost.mirrors));
    end
end

%% plot: velocity potential
Z_copy = Z/max(max(Z));
figure(fig1);
clf
imshow(Z_copy); % show the "image"
colormap spring;

%% plot: ghosts and mirror points
K = zeros(size(Z));
K(vertcat(ghosts_struct.mirrors)) = 1;
K([vertcat(ghosts_struct.row),vertcat(ghosts_struct.col)]) = 0.5;
figure(fig1);
clf
imshow(K);

%% plot: animation and velocity vector field streamlines
[dx, dy] = gradient(Z);

% discard vectors inside or at the boundary of the shape
dx(interior | ghosts) = 0;
dy(interior | ghosts) = 0;

% crop out vectors near the discontinuities
cropx = [-2 6];
cropy = [-3 3];
crop_mat = cropx(1)<=xx & xx<cropx(2) & cropy(1)<yy & yy<cropy(2);
crop = @(x) x(any(crop_mat,2),any(crop_mat,1));
dx = -crop(dx);
dy = -crop(dy);
xxcrop = crop(xx);
yycrop = crop(yy);
Zcrop = crop(Z);

% make Z occupy full color range
Zcrop = Zcrop*(255/max(max(Zcrop))); 

% configure figure
figure(fig2)
clf
daspect([1 1 1])

% background 
image('CData',Zcrop,'XData',cropx,'YData',cropy);
colormap spring;
axis equal

% particle animation
[verts,averts] = streamslice(xxcrop,yycrop,dx,dy); 
sl = streamline([verts averts]);
axis tight manual off equal;
ax = gca;
ax.Position = [0,0,1,1];
set(sl,'Visible','on')
iverts = interpstreamspeed(xxcrop,yycrop,dx,dy,verts,.5);
streamparticles(iverts, 200, ...
    'Animate',205,'FrameRate',40, ...
    'MarkerSize',10,'MarkerFaceColor',[1 1 1])

%% plot: velocity field streamlines only
% streamplot starts at x=0, distributed along y
stream_spacing = 0.4;
startys=cropy(1):stream_spacing:cropy(2);
startxs=ones(size(startys))*cropx(1);

figure(fig2)
clf
quiver(xxcrop,yycrop,dx,dy);
streamline(xxcrop,yycrop,dx,dy,startxs,startys);
axis equal

%% functions
function [exterior, ghosts, ghosts_struct, interior] = elliptic_curve(xmin,xmax,ymin,ymax,imagepath,x,y)
    [xx,yy] = meshgrid(x,y);
    
    shape = imread(imagepath);
    shape = shape(:,:,1) < 123; % all black pixels
    width = xmax-xmin;
    height = ymax-ymin;
    shp_dx = width / size(shape,2);
    shp_dy = height / size(shape,1);
    
    Z = false(size(xx));
    for row = 1:size(shape,1)
        for col = 1:size(shape,2)
            if shape(row,col) %inside the shape
                shp_x = xmin + shp_dx*col;
                shp_y = ymin + shp_dy*row;
                dists = sqrt(abs(xx-shp_x).^2 + abs(yy-shp_y).^2);
                [~,I] = min(dists(:));
                Z(I) = 1;
            end
        end
    end
    
    Z_n = north(Z);
    Z_e = east(Z);
    Z_s = south(Z);
    Z_w = west(Z);
    
    interior = Z_n & Z_e & Z_s & Z_w;
    exterior = ~Z;
    
    ghosts = xor(Z,interior);
    ghosts_struct = [];
    [ghost_rows, ghost_cols] = find(ghosts);
    for i = length(ghost_rows):-1:1 %backwards to preallocate
        ix = ghost_rows(i);
        iy = ghost_cols(i);
        stenciled = false(size(Z));
        stenciled(stencil_25(iy,ix,Z)) = 1;
        stenciled = stenciled & exterior;
        stenciled = find(stenciled);
        ghosts_struct(i).row = ix;
        ghosts_struct(i).col = iy;
        ghosts_struct(i).mirrors = stenciled;
    end
end

% outputs indicies of 25-point stencil around point (i,j)
function out = stencil_25(i,j,grid)
    is = i-2:i+2;
    js = j-2:j+2;
    [ii, jj] = meshgrid(is,js);
    out = reshape((ii-1)*size(grid,1)+jj,[],1);
end

% outputs indicies of 9-point stencil around point (i,j)
function out = stencil_9(i,j,grid)
    is = i-1:i+1;
    js = j-1:j+1;
    [ii, jj] = meshgrid(is,js);
    out = reshape((ii-1)*size(grid,1)+jj,[],1);
end

function out = north(a)
    out = padarray(a(1:end-1,:),[1,0],'pre');
end

function out = east(a)
    out = padarray(a(:,2:end),[0,1],'post');
end

function out = south(a)
    out = padarray(a(2:end,:),[1,0],'post');
end

function out = west(a)
    out = padarray(a(:,1:end-1),[0,1],'pre');
end
