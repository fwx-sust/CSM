function [val_xi,val_eta] = Tri_grad(spl,aa)
if spl==1
    if aa == 1
        val_xi  = -0.5;
        val_eta = -0.5;
    elseif aa == 2
        val_xi=0.5;
        val_eta=0;
    elseif aa == 3
        val_xi=0;
        val_eta=0.5;
    else
        error('Error: value of a should be 1,2,3');
    end
end

if spl==2
    if aa == 1
        val_xi=0;
        val_eta=-0.5;
    elseif aa == 2
        val_xi=0.5;
        val_eta=0.5;
    elseif aa == 3
        val_xi=-0.5;
        val_eta=0;
    else
        error('Error: value of a should be 1,2,3');
    end
end