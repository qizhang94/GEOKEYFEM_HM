function snscontour(node, element,fac, u_xp,u_yp,variable,name)
    node2=node+fac*[u_xp u_yp];
    for k=1:length(element(:,1))
        i=element(k,:);  % Connectivity Matrix
        x=node2(i,1);
        y=node2(i,2);
        s=variable(i);
        h2=fill(x,y,s,'FaceColor','interp'); hold on
    end
    %%
    shading interp;  % mesh or not
    % shading flat
    colormap(jet)
    % colormap(cmapRainbow);
    c = colorbar;
    %%
    c.Label.String =name;
    brighten(0.3);
    %%
    axis off
    box off
    axis equal
end