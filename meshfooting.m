function [node,element]=meshfooting(mesh_type, keypoint, meshsize, inner_edge)
addpath(genpath('./dengwirda-mesh2d-ceb68eb'))
addpath(genpath('./mesh2D_MLC'))
addpath('./meshingN')
addpath('./plottingN')
warning off
%%
global geoheight
global geowidth
global hhmax
global hhmin
global refine_x
global refine_y

geoheight = keypoint(3,2);
geowidth = keypoint(3,1);
if inner_edge
    edge0 = [1,2; 2,3; 3,4; 4,1; 5,6; 6,7];
else
    edge0 = [1,2; 2,3; 3,4; 4,1];
end

part{1} = [1,2,3,4];
%---------------------------------------------- do size-fun.
hhmax = meshsize(1);
hhmin = meshsize(2);
if length(meshsize) == 2
    refine_x = 0.1*geowidth;
    refine_y = 0.7*geoheight;
else
    refine_x = meshsize(3);
    refine_y = meshsize(4);
end

hfun = @hfun_footing;
%---------------------------------------------- do mesh-gen.
opts.kind = 'delaunay';
switch mesh_type
    case 'Uniform mesh'
        [vert, etri,tria, tnum] = refine2(keypoint, edge0, part, opts, hhmax) ;
    case 'Non-uniform mesh'
        [vert, etri,tria, tnum] = refine2(keypoint, edge0, part, opts, hfun);
end
[node, ~, element, ~] = smooth2(vert,etri,tria,tnum);
end


function h = hfun_footing(node)
global geoheight
global geowidth
global hhmax
global hhmin
global refine_x
global refine_y
%% User defined size function for square
Hmax=hhmax;
Hmin=hhmin;
%
x=node(:, 1);
y=node(:, 2);
%
A=geowidth;   % max(x)-min(x); length in x 
B=geoheight;   % max(y)-min(y); length in y 
%
C=((x-refine_x)+abs(refine_x-x))/2;
%
D=((refine_y-y)+abs(y-refine_y))/2;
%
h=(Hmax-Hmin).*((C/A).^2+(D/B).^2)+Hmin;
end