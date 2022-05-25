function residual_traction = assign_tractionBC2(residual_traction, b_top, traction_f, time)

global gcoord
ndof = 2;
residual_traction(:) = 0;

temp = [gcoord(b_top, 1), b_top];
temp = sortrows(temp);
b_top = temp(:, 2);


for i = 1:length(b_top)
    if i==1
        a = gcoord(b_top(1), 1);
        b = gcoord(b_top(2), 1);
        fun = @(x)((b-x)/(b-a).*traction_f(x, time));
        residual_traction(ndof*b_top(1)) = integral(fun, a, b);

    elseif i == length(b_top)
        a = gcoord(b_top(end-1), 1);
        b = gcoord(b_top(end), 1);
        fun = @(x)((x-a)/(b-a).*traction_f(x, time));
        residual_traction(ndof*b_top(end)) = integral(fun, a, b);

    else
        a = gcoord(b_top(i-1), 1);
        b = gcoord(b_top(i+1), 1);
        c = gcoord(b_top(i), 1);
        fun = @(x)((x-a)/(c-a).*traction_f(x, time));
        residual_traction(ndof*b_top(i)) = residual_traction(ndof*b_top(i)) + integral(fun, a, c);
        fun = @(x)((b-x)/(b-c).*traction_f(x, time));
        residual_traction(ndof*b_top(i)) = residual_traction(ndof*b_top(i)) + integral(fun, c, b);

    end
end

end