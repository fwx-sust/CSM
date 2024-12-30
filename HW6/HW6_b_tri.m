%triangle element
clear all; clc;

kappa = 1.0; % conductivity

% exact solution
exact = @(x,y) x*(1-x)*y*(1-y);
exact_x = @(x,y) (1-2*x)*y*(1-y);
exact_y = @(x,y) x*(1-x)*(1-2*y);

f = @(x,y) 2.0*kappa*x*(1-x) + 2.0*kappa*y*(1-y); % source term

% quadrature rule
n_int_xi  = 3;
n_int_eta = 3;
n_int     = n_int_xi * n_int_eta;
[xi, eta, weight] = Gauss2D(n_int_xi, n_int_eta);

% mesh generation
n_en   = 3;               % number of nodes in an element
n_el_x0=[20,30,40,50,60,70,80,90]*2; %changable
n_el_y0=[20,30,40,50,60,70,80,90]*2; %changable
flag=0;
for n_el_x = n_el_x0              % number of elements in x-dir
    for n_el_y = n_el_y0               % number of elements in y-dir
        n_el   = 0.5*n_el_x * n_el_y; % total number of elements: triangles

        n_np_x = 0.5*n_el_x + 1;      % number of nodal points in x-dir
        n_np_y = 0.5*n_el_y + 1;      % number of nodal points in y-dir
        n_np   = n_np_x * n_np_y; % total number of nodal points

        x_coor = zeros(n_np, 1);
        y_coor = x_coor;

        hx = 2.0 / n_el_x;        % mesh size in x-dir
        hy = 2.0 / n_el_y;        % mesh size in y-dir
        h0=sqrt(hx^2 +hy^2);
        if flag==0
            flag=flag+1;
            hh(flag)=h0;
        else
            %log(e)-log(h) is a function
            if ismember(h0,hh)==1
                continue;
            else
                flag=flag+1;
                hh(flag)=h0;
            end
        end

        % generate the nodal coordinates
        for ny = 1 : n_np_y
            for nx = 1 : n_np_x
                index = (ny-1)*n_np_x + nx; % nodal index
                x_coor(index) = (nx-1) * hx;
                y_coor(index) = (ny-1) * hy;
            end
        end

        % IEN array
        IEN = zeros(n_el, n_en);
        for ex = 1 : 0.5*n_el_x
            for ey = 1 : 0.5*n_el_y
                for sp=1:2 % a quadrilateral is subdivided into two triangle
                    ee = 2*( (ey-1) * 0.5*n_el_x + ex )- (2-sp); % element index
                    if sp==1
                        IEN(ee, 1) = (ey-1) * n_np_x + ex;
                        IEN(ee, 2) = (ey-1) * n_np_x + ex + 1;
                        IEN(ee, 3) =  ey    * n_np_x + ex + 1;
                    end
                    if sp==2
                        IEN(ee, 1) = (ey-1) * n_np_x + ex+1;
                        IEN(ee, 2) =  ey    * n_np_x + ex + 1;
                        IEN(ee, 3) =  ey    * n_np_x + ex;
                    end
                end
            end
        end

        % ID array
        ID = zeros(n_np,1);
        counter = 0;
        for ny = 2 : n_np_y - 1
            for nx = 2 : n_np_x - 1
                index = (ny-1)*n_np_x + nx;
                counter = counter + 1;
                ID(index) = counter;
            end
        end
        n_eq = counter;
        LM = ID(IEN);

        % allocate the stiffness matrix and load vector
        K = spalloc(n_eq, n_eq, 9 * n_eq);
        F = zeros(n_eq, 1);

        % loop over element to assembly the matrix and vector
        for ee = 1 : n_el
            if mod(ee,2)==0
                sp=2;
            else
                sp=1;
            end
            x_ele = x_coor( IEN(ee, 1:n_en) );
            y_ele = y_coor( IEN(ee, 1:n_en) );

            k_ele = zeros(n_en, n_en); % element stiffness matrix
            f_ele = zeros(n_en, 1);    % element load vector

            for ll = 1 : n_int
                x_l = 0.0; y_l = 0.0;
                dx_dxi = 0.0; dx_deta = 0.0;
                dy_dxi = 0.0; dy_deta = 0.0;
                for aa = 1 : n_en
                    x_l = x_l + x_ele(aa) * Tri(sp,aa, xi(ll), eta(ll));
                    y_l = y_l + y_ele(aa) * Tri(sp,aa, xi(ll), eta(ll));
                    [Na_xi, Na_eta] = Tri_grad(sp,aa);
                    dx_dxi  = dx_dxi  + x_ele(aa) * Na_xi;
                    dx_deta = dx_deta + x_ele(aa) * Na_eta;
                    dy_dxi  = dy_dxi  + y_ele(aa) * Na_xi;
                    dy_deta = dy_deta + y_ele(aa) * Na_eta;
                end

                detJ = dx_dxi * dy_deta - dx_deta * dy_dxi;

                for aa = 1 : n_en
                    Na = Tri(sp,aa, xi(ll), eta(ll));
                    [Na_xi, Na_eta] = Tri_grad(sp,aa);
                    Na_x = (Na_xi * dy_deta - Na_eta * dy_dxi) / detJ;
                    Na_y = (-Na_xi * dx_deta + Na_eta * dx_dxi) / detJ;

                    f_ele(aa) = f_ele(aa) + weight(ll) * detJ * f(x_l, y_l) * Na;

                    for bb = 1 : n_en
                        Nb = Tri(sp,bb, xi(ll), eta(ll));
                        [Nb_xi, Nb_eta] = Tri_grad(sp,bb);
                        Nb_x = (Nb_xi * dy_deta - Nb_eta * dy_dxi) / detJ;
                        Nb_y = (-Nb_xi * dx_deta + Nb_eta * dx_dxi) / detJ;

                        k_ele(aa, bb) = k_ele(aa,bb) + weight(ll) * detJ * kappa * (Na_x * Nb_x + Na_y * Nb_y);
                    end % end of bb loop
                end % end of aa loop
            end % end of quadrature loop

            for aa = 1 : n_en
                PP = LM(ee, aa);
                if PP > 0
                    F(PP) = F(PP) + f_ele(aa);

                    for bb = 1 : n_en
                        QQ = LM(ee, bb);
                        if QQ > 0
                            K(PP, QQ) = K(PP, QQ) + k_ele(aa, bb);
                        else
                            % modify F with the boundary data
                            % here we do nothing because the boundary data g is zero or
                            % homogeneous
                        end
                    end
                end
            end
        end

        % solve the stiffness matrix
        dn = K \ F;

        % insert dn back into the vector for all nodes
        disp = zeros(n_np, 1);

        for ii = 1 : n_np
            index = ID(ii);
            if index > 0
                disp(ii) = dn(index);
            else
                % modify disp with the g data. Here it does nothing because g is zero
            end
        end

        fenzi0=0;
        fenzi1=0;
        for ee=1:n_el
            if mod(ee,2)==0
                sp=2;
            else
                sp=1;
            end
            x_ele = x_coor( IEN(ee, 1:n_en) );
            y_ele = y_coor( IEN(ee, 1:n_en) );
            u_ele = disp( IEN(ee, :) );
            for qua=1:n_int
                %calculate error0 and error1
                uh=0;u=0;
                x=0;y=0;
                uhx=0;ux=0;uhy=0;
                dx_dxi = 0.0; dx_deta = 0.0;
                dy_dxi = 0.0; dy_deta = 0.0;

                for aa = 1 : n_en
                    uh=uh+u_ele(aa)*Tri(sp,aa, xi(qua), eta(qua)); %qN
                    x=x+x_ele(aa)*Tri(sp,aa, xi(qua), eta(qua));
                    y=y+y_ele(aa)*Tri(sp,aa, xi(qua), eta(qua));
                    u=exact(x,y);

                    [Na_xi, Na_eta] = Tri_grad(sp,aa);
                    %because it is relative to sp, so J is different.
                    dx_dxi  = dx_dxi  + x_ele(aa) * Na_xi;
                    dx_deta = dx_deta + x_ele(aa) * Na_eta;
                    dy_dxi  = dy_dxi  + y_ele(aa) * Na_xi;
                    dy_deta = dy_deta + y_ele(aa) * Na_eta;

                end
                detJ = dx_dxi * dy_deta - dx_deta * dy_dxi;

                for aa= 1: n_en
                    [Na_xi, Na_eta] = Tri_grad(sp,aa);
                    Na_x = (Na_xi * dy_deta - Na_eta * dy_dxi) / detJ;
                    Na_y = (-Na_xi * dx_deta + Na_eta * dx_dxi) / detJ;
                    uhx=uhx+u_ele(aa)*(Na_xi * dy_deta - Na_eta * dy_dxi) / detJ;
                    uhy=uhy+u_ele(aa)*(-Na_xi * dx_deta + Na_eta * dx_dxi) / detJ ;
                    ux=exact_x(x,y);
                    uy=exact_y(x,y);
                end
                %page174: c=0.5 for triangles
                fenzi0= fenzi0+ 0.5*detJ*weight(qua)*(uh-u)^2;
                fenzi1= fenzi1+ 0.5*detJ*weight(qua)*( (uh-u)^2 + (uhx -ux)^2 + (uhy -uy)^2 );
            end
        end
        fenzi0=sqrt(fenzi0);
        error0(flag)=fenzi0;
        fenzi1=sqrt(fenzi1);
        error1(flag)=fenzi1;
    end
end
%plot
figure()
plot(log(hh),log(error0),'ro');
hold on
A(1,:)=polyfit(log(hh),log(error0),1);
xa=linspace(-5,-2,20);ya=A(1,1)*xa+A(1,2);
plot(xa,ya,'-');

plot(log(hh),log(error1),'x');
hold on
B(1,:)=polyfit(log(hh),log(error1),1);
xb=linspace(-5,-2,20);yb=B(1,1)*xb+B(1,2);
plot(xb,yb,'-')
