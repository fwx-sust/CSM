%error calculator————可以用来观察网格大小固定的情况。
fenzi0_x=0;
fenzi1_x=0;
fenzi0_y=0;
fenzi1_y=0;
for ee=1:n_el
    x_ele = x_coor( IEN(ee, 1:n_en) );
    y_ele = y_coor( IEN(ee, 1:n_en) );
    ux_ele = disp( IEN(ee, :),1 );
    uy_ele = disp( IEN(ee, :),2 );
    for qua=1:n_int
        %calculate error0 and error1
        x=0;y=0;
        uhx=0;uhy=0;
        uhx_x=0;uhx_y=0;uhy_x=0;uhy_y=0;
        ux=0;uy=0;
        ux_x=0;ux_y=0;uy_x=0;uy_y=0;

        for aa = 1 : n_en
            uhx=uhx+ux_ele(aa)*Quad(aa, xi(qua), eta(qua));
            uhy=uhy+uy_ele(aa)*Quad(aa, xi(qua), eta(qua));
            x=x+x_ele(aa)*Quad(aa, xi(qua), eta(qua));
            y=y+y_ele(aa)*Quad(aa,xi(qua),eta(qua));
            ux=fux(x,y);
            uy=fuy(x,y);

            [Na_xi, Na_eta] = Quad_grad(aa, xi(qua), eta(qua));
            uhx_x=uhx_x+ux_ele(aa)*(Na_xi * dy_deta - Na_eta * dy_dxi) / detJ;
            uhx_y=uhx_y+ux_ele(aa)*(-Na_xi * dx_deta + Na_eta * dx_dxi) / detJ ;
            uhy_x=uhy_x+uy_ele(aa)*(Na_xi * dy_deta - Na_eta * dy_dxi) / detJ;
            uhy_y=uhy_y+uy_ele(aa)*(-Na_xi * dx_deta + Na_eta * dx_dxi) / detJ ;
            ux_x=fux_x(x,y);
            ux_y=fux_y(x,y);
            uy_x=fuy_x(x,y);
            uy_y=fuy_y(x,y);
        end
        fenzi0_x= fenzi0_x+detJ*weight(qua)*(uhx-ux)^2;
        fenzi0_y= fenzi0_y+detJ*weight(qua)*(uhy-uy)^2;

        fenzi1_x= fenzi1_x+ detJ*weight(qua)*( (uhx-ux)^2 + (uhx_x -ux_x)^2 + (uhx_y -ux_y)^2 );
        fenzi1_y= fenzi1_y+ detJ*weight(qua)*( (uhy-uy)^2 + (uhy_x -uy_x)^2 + (uhy_y -uy_y)^2 );
    end
end
fenzi0_x=sqrt(fenzi0_x);
error0_x(flag)=fenzi0_x;
fenzi0_y=sqrt(fenzi0_y);
error0_y(flag)=fenzi0_y;

fenzi1_x=sqrt(fenzi1_x);
error1_x(flag)=fenzi1_x;
fenzi1_y=sqrt(fenzi1_y);
error1_y(flag)=fenzi1_y;