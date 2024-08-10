function t = smooth_result(node, field, num_gridx, num_gridy, filter)
%SMOOTH_RESULT 光滑模拟结果云图
%   通过散点插值尽可能消除结果的锯齿或者震荡
%   https://www.mathworks.com/help/matlab/ref/griddata.html#mw_0ad56d6b-0679-4a71-a305-8196feb4fd7a
%   [Xq,Yq,vq] = griddata(x,y,v,xq,yq,method)
%   前提是模拟的区域是一个规则的长方形 reference configuration

field(field<filter(1)) = filter(1);
field(field>filter(2)) = filter(2);

[xq,yq] = meshgrid(...
            linspace(min(node(:,1)),max(node(:,1)),num_gridx),...
            linspace(min(node(:,2)),max(node(:,2)),num_gridy)...
          );
F = scatteredInterpolant(node(:,1),node(:,2),field,'natural','nearest');
vq = F(xq, yq);

%[Xq,Yq,vq] = griddata(node(:,1),node(:,2),field,xq,yq,method); % "natural" or "v4", do not use "cubic"!

[~,h]=contourf(xq,yq,vq,100);
set(h, 'edgecolor','none');
colormap("jet");

% AVAQUS_CM_HTML = {'#FF0000', '#FF5D00', '#FFB900','#E8FF00', '#8BFF00',...
%     '#2EFF00', '#00FF2E', '#00FF8B', '#00FFE8', '#00B9EF', '#005DEF', '#0000FF'};
% mycolormap = customcolormap(linspace(0,1,12), AVAQUS_CM_HTML, 32);
% colormap(mycolormap);

c = colorbar;
t = get(c,'Limits');
axis equal;
end