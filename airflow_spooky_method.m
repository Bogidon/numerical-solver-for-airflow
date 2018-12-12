%% Lin II Final: Air Flow Around Elliptical Objects
% old method – kept for archiving, do not use

%% figures
fig1 = figure;
fig2 = figure;

%% the code

% create the scalar grid
spacing = 0.1;
x = -10:spacing:10;
y = -20:spacing:20;
[xx,yy] = meshgrid(x,y);
[exterior, surface, ghosts, mirrors, interior] = elliptic_curve(-3,0,3,2,x,y);
left_side = logical(padarray(ones(length(y),1),[0,length(x)-1],'post'));
right_side = logical(padarray(ones(length(y),1),[0,length(x)-1],'pre'));
top_side = logical(padarray(ones(1,length(x)),[length(y)-1,0],'post'));
bottom_side = logical(padarray(ones(1,length(x)),[length(y)-1,0],'pre'));
sides = left_side | right_side |top_side | bottom_side;
exterior = xor(exterior,sides);
Z = zeros(size(exterior)) + (left_side * 20);

for i = 1:1000000
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
    for k = 1:length(ghosts)
        ghost = ghosts(k);
        mirror = mirrors(k);
        Z(ghost) = Z(mirror);
    end
end
%%
Z = Z/4;
figure(fig1);
clf
imshow(Z); % show the "image"
colormap spring;
%% gradient stuff
crop_mat = -7<xx & xx<3 & -4<yy & yy<4;
crop = @(x) x(any(crop_mat,2),any(crop_mat,1));
[dx, dy] = gradient(Z);
dx = -crop(dx);
dy = -crop(dy);
xxcrop = crop(xx);
yycrop = crop(yy);

figure(fig2)
clf
quiver(xxcrop,yycrop,dx,dy)

%% functions
function [exterior, surface, ghosts, mirrors, interior] = elliptic_curve(h,k,a,b,x,y)
    [xx,yy] = meshgrid(x,y);
    Z=(((xx-h).^2)/a^2)+(((yy-k).^2)/b^2);
    Z=Z<=1;
    
    Z_n = north(Z);
    Z_ne = north(east(Z));
    Z_e = east(Z);
    Z_se = south(east(Z));
    Z_s = south(Z);
    Z_sw = south(west(Z));
    Z_w = west(Z);
    Z_nw = north(west(Z));
    
    interior = Z_n & Z_e & Z_s & Z_w;
    surface = xor(Z,interior);
    
    S_n = north(surface);
    S_ne = north(east(surface));
    S_e = east(surface);
    S_se = south(east(surface));
    S_s = south(surface);
    S_sw = south(west(surface));
    S_w = west(surface);
    S_nw = north(west(surface));
    
    spooky_surface_n = (S_e & S_w & Z_s);
    spooky_surface_ne = (S_nw & S_se & Z_sw);
    spooky_surface_e = (S_n & S_s & Z_w);
    spooky_surface_se = (S_ne & S_sw & Z_nw);
    spooky_surface_s = (S_e & S_w & Z_n);
    spooky_surface_sw = (S_nw & S_se & Z_ne);
    spooky_surface_w = (S_n & S_s & Z_e);
    spooky_surface_nw = (S_ne & S_sw & Z_se);
    
    ghosts = find(...
        south(spooky_surface_n)...
        | south(west(spooky_surface_ne))...
        | west(spooky_surface_e)...
        | north(west(spooky_surface_se))...
        | north(spooky_surface_s)...
        | north(east(spooky_surface_sw))...
        | east(spooky_surface_w)...
        | south(east(spooky_surface_nw)));
    
    mirrors = find(...
        north(spooky_surface_n)...
        | north(east(spooky_surface_ne))...
        | east(spooky_surface_e)...
        | south(east(spooky_surface_se))...
        | south(spooky_surface_s)...
        | south(west(spooky_surface_sw))...
        | west(spooky_surface_w)...
        | north(west(spooky_surface_nw)));
    
%     ghosts_and_mirrors = [find(ghosts) find(mirrors)];
%     ghosts_and_mirrors = ghosts | mirrors;
    exterior = ~Z;
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
