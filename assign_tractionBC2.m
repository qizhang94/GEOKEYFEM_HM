function residual_traction = assign_tractionBC2(residual_traction, b_natural, traction_f, time, boundary)
global gcoord
ndof = 2;
residual_traction(:) = 0;

if nargin <= 4
    boundary = 'horizontal';
end

if strcmp(boundary, 'horizontal')
    temp = [gcoord(b_natural, 1), b_natural];
    temp = sortrows(temp);
    b_natural = temp(:, 2);
    
    for i = 1:length(b_natural)
        if i==1
            a = gcoord(b_natural(1), 1);
            b = gcoord(b_natural(2), 1);
            fun = @(x)((b-x)/(b-a).*traction_f(x, time));
            residual_traction(ndof*b_natural(1)-1:ndof*b_natural(1)) = integral(fun, a, b,'ArrayValued',true);
    
        elseif i == length(b_natural)
            a = gcoord(b_natural(end-1), 1);
            b = gcoord(b_natural(end), 1);
            fun = @(x)((x-a)/(b-a).*traction_f(x, time));
            residual_traction(ndof*b_natural(end)-1:ndof*b_natural(end)) = integral(fun, a, b, 'ArrayValued',true);
    
        else
            a = gcoord(b_natural(i-1), 1);
            b = gcoord(b_natural(i+1), 1);
            c = gcoord(b_natural(i), 1);
            fun = @(x)((x-a)/(c-a).*traction_f(x, time));
            residual_traction(ndof*b_natural(i)-1:ndof*b_natural(i)) = ...
                residual_traction(ndof*b_natural(i)-1:ndof*b_natural(i)) + integral(fun, a, c, 'ArrayValued',true);
            fun = @(x)((b-x)/(b-c).*traction_f(x, time));
            residual_traction(ndof*b_natural(i)-1:ndof*b_natural(i)) = ...
                residual_traction(ndof*b_natural(i)-1:ndof*b_natural(i)) + integral(fun, c, b, 'ArrayValued',true);
    
        end
    end
elseif strcmp(boundary, 'vertical')
    temp = [gcoord(b_natural, 2), b_natural]; % 2 columns
    temp = sortrows(temp);
    b_natural = temp(:, 2);
    
    for i = 1:length(b_natural)
        if i==1
            a = gcoord(b_natural(1), 2);
            b = gcoord(b_natural(2), 2);
            fun = @(x)((b-x)/(b-a).*traction_f(x, time));
            residual_traction(ndof*b_natural(1)-1:ndof*b_natural(1)) = integral(fun, a, b, 'ArrayValued',true);
    
        elseif i == length(b_natural)
            a = gcoord(b_natural(end-1), 2);
            b = gcoord(b_natural(end), 2);
            fun = @(x)((x-a)/(b-a).*traction_f(x, time));
            residual_traction(ndof*b_natural(end)-1:ndof*b_natural(end)) = integral(fun, a, b, 'ArrayValued',true);
    
        else
            a = gcoord(b_natural(i-1), 2);
            b = gcoord(b_natural(i+1), 2);
            c = gcoord(b_natural(i), 2);
            fun = @(x)((x-a)/(c-a).*traction_f(x, time));
            residual_traction(ndof*b_natural(i)-1:ndof*b_natural(i)) = ...
                residual_traction(ndof*b_natural(i)-1:ndof*b_natural(i)) + integral(fun, a, c, 'ArrayValued',true);
            fun = @(x)((b-x)/(b-c).*traction_f(x, time));
            residual_traction(ndof*b_natural(i)-1:ndof*b_natural(i)) = ...
                residual_traction(ndof*b_natural(i)-1:ndof*b_natural(i)) + integral(fun, c, b, 'ArrayValued',true);
    
        end
    end
end

end