function u = u_ol(t, h, u0)
    if rem(t, h) < .00001
        % We're at a time which is a whole multiple of h
        u = u0(round(t/h));
    elseif rem(t*2, h) < .00001
        % We're at a midpoint (knot point)
        n = round((t-h)/h);
        u = (u0(n)+u0(n+1))/2;
    else
        error('Not a whole or half multiple of h');
    end
end