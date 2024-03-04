function cstate = astate2cstate(astate, xgrid, ygrid)
    w = size(xgrid, 2);
    i = ceil(astate / w);
    j = astate - (i - 1)  * w;
    cstate = [xgrid(i,j); ygrid(i,j)];  
end

