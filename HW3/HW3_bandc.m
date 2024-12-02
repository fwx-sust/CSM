%question b and c
%Calculate the following two relative errors of the solution using your code

clear; clc; clf; % clean the memory, screen, and figure

% Problem definition
f = @(x) -20*x.^3; % f(x) is the source
g = 1.0;           % u    = g  at x = 1
h = 0.0;           % -u,x = h  at x = 0

% Setup the mesh
pd   = [1,2,3];              % polynomial degree
for pp=pd
    n_en = pp + 1;         % number of element or local nodes
    nel=[2,4,6,8,10,12,14,16];
    error_L2 = zeros(size(nel));
    error_H1 = zeros(size(nel));
    n_int = 10;
    for n_el = nel              %2, 4, 6, 8, 10, 12, 14, and 16 elements
        n_np = n_el * pp + 1;  % number of nodal points
        n_eq = n_np - 1;       % number of equations
        fenziL2=0;
        fenmuL2=0;
        fenziH1=0;
        fenmuH1=0;

        hh(n_el/2) = 1.0 / (n_np - 1); % space between two adjacent nodes
        x_coor = 0 : hh(n_el/2) : 1;   % nodal coordinates for equally spaced nodes

        IEN = zeros(n_el, n_en);

        for ee = 1 : n_el
            for aa = 1 : n_en
                IEN(ee, aa) = (ee - 1) * pp + aa;
            end
        end

        % Setup the ID array for the problem
        ID = 1 : n_np;
        ID(end) = 0;

        % Setup the quadrature rule
        [xi, weight] = Gauss(n_int, -1, 1);

        % allocate the stiffness matrix
        K = spalloc(n_eq, n_eq, (2*pp+1)*n_eq);
        F = zeros(n_eq, 1);

        % Assembly of the stiffness matrix and load vector
        for ee = 1 : n_el
            k_ele = zeros(n_en, n_en); % allocate a zero element stiffness matrix
            f_ele = zeros(n_en, 1);    % allocate a zero element load vector

            x_ele = x_coor(IEN(ee,:)); % x_ele(aa) = x_coor(A) with A = IEN(aa, ee)

            % quadrature loop
            for qua = 1 : n_int
                dx_dxi = 0.0; %differential
                x_l = 0.0;

                for aa = 1 : n_en
                    x_l    = x_l    + x_ele(aa) * PolyShape(pp, aa, xi(qua), 0);
                    dx_dxi = dx_dxi + x_ele(aa) * PolyShape(pp, aa, xi(qua), 1);
                end

                dxi_dx = 1.0 / dx_dxi;

                for aa = 1 : n_en
                    f_ele(aa) = f_ele(aa) + weight(qua) * PolyShape(pp, aa, xi(qua), 0) * f(x_l) * dx_dxi;
                    for bb = 1 : n_en
                        k_ele(aa, bb) = k_ele(aa, bb) + weight(qua) * PolyShape(pp, aa, xi(qua), 1) * PolyShape(pp, bb, xi(qua), 1) * dxi_dx;
                    end
                end
            end

            % Assembly of the matrix and vector based on the ID or LM data
            for aa = 1 : n_en
                P = ID(IEN(ee,aa));
                if(P > 0)
                    F(P) = F(P) + f_ele(aa);
                    for bb = 1 : n_en
                        Q = ID(IEN(ee,bb));
                        if(Q > 0)
                            K(P, Q) = K(P, Q) + k_ele(aa, bb);
                        else
                            F(P) = F(P) - k_ele(aa, bb) * g; % handles the Dirichlet boundary data
                        end
                    end
                end
            end
        end

        % ee = 1 F = NA(0)xh
        F(ID(IEN(1,1))) = F(ID(IEN(1,1))) + h;

        % Solve Kd = F equation
        d_temp = K \ F;
        disp = [d_temp; g];

        for ee = 1 : n_el
            x_ele = x_coor( IEN(ee, :) );
            u_ele = disp( IEN(ee, :) );

            if ee == n_el
                n_sam_big = 2;
            else
                n_sam_big =1;
            end

            for i = n_sam_big : n_int
                for qua=1:n_int
                    %calculate errorL2
                    uh=0;x=0;u=0;
                    for aa = 1 : n_en
                        uh=uh+u_ele(aa)*PolyShape(pp, aa, xi(i), 0);
                        x=x+x_ele(aa)*PolyShape(pp, aa, xi(i), 0);
                        u=x.^5;
                    end
                    %一个点
                    fenziL2= fenziL2+ weight(qua)*(uh-u)^2;
                    fenmuL2=fenmuL2+weight(qua)* (u^2);

                    %calculate errorH1
                    uhx=0;xx=0;ux=0;
                    for aa = 1 : n_en
                        uhx=uhx+u_ele(aa)*PolyShape(pp, aa, xi(i), 1);
                        xx=xx+x_ele(aa)*PolyShape(pp, aa, xi(i), 0);
                        ux=5* (xx.^4);
                    end
                    fenziH1= fenziH1+ weight(qua)*(uhx-ux)^2;
                    fenmuH1=fenmuH1+weight(qua)* (ux^2);

                end
            end
        end
        %21个点的值相加
        fenziL2=sqrt(fenziL2);
        fenmuL2=sqrt(fenmuL2);
        error_L2(n_el/2) = fenziL2/fenmuL2;
        fenziH1=sqrt(fenziH1);
        fenmuH1=sqrt(fenmuH1);
        error_H1(n_el/2) = fenziH1/fenmuH1;
    end
    figure()
    plot(log(hh),log(error_L2),'ro');
    hold on
    title(pp,"degree's L2")
    A(pp,:)=polyfit(log(hh),log(error_L2),1);
    xa=linspace(-4,0,20);ya=A(pp,1).*xa+A(pp,2);
    plot(xa,ya,'-k')
    figure()
    plot(log(hh),log(error_H1),'x');
    title(pp,"degree's H1")
    hold on
    B(pp,:)=polyfit(log(hh),log(error_H1),1);
    xb=linspace(-4,0,20);yb=B(pp,1).*xb+B(pp,2);
    plot(xb,yb,'-k')
end