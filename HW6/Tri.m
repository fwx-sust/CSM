%Triangle
function val = Tri(spl,aa, xi, eta)
if spl==1
    if aa == 1
        val = -0.5*xi-0.5*eta;
    elseif aa == 2
        val = 0.5+0.5*xi;
    elseif aa == 3
        val = 0.5+0.5*eta;
    else
        error('Error: value of a should be 1,2,3');
    end
end

if spl==2
    if aa == 1
        val = 0.5-0.5*eta;
    elseif aa == 2
        val = 0.5*(xi+eta);
    elseif aa == 3
        val = 0.5-0.5*xi;
    else
        error('Error: value of a should be 1,2,3');
    end
end