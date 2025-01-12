%error for 1*1矩形
%L2
%quadrilateral element
clear all; clc;

E=2e11; %Young's modulus
v=0.3; %Poisson's ratio

%D(x)————specify either plane strain or plane stress
question_def = 1; %1 for plane strain, 2 for plane stress
D0 = stiffnessD(question_def, E, v); %stress=D*strain

%body force
%只适用于plane strain，因为应变和应力的关系为plane strain提供的
Da=E/(v+1)/(2*v-1);
f=@(x,y,dof) (dof==1)*(-Da)*( 2*(v-1)*(y^2-y)-v*(2*x-1)*(2*y-1)+(2*v-1)*(x^2-x)+0.5*(2*v-1)*(2*x-1)*(2*y-1) ) +...
    (dof==2)*(-Da)*( 2*(v-1)*(x^2-x)-v*(2*y-1)*(2*x-1)+(2*v-1)*(y^2-y)+0.5*(2*v-1)*(2*y-1)*(2*x-1) ) ;

%exact solution————以下是最简单的两端固支的位移
%位移理论值
fux=@(x,y) x.*(x-1).*y.*(y-1); %点乘方便后续画图
fuy=@(x,y) x.*(x-1).*y.*(y-1);

fux_x=@(x,y) (2.*x-1).*(y.^2-y);
fux_y=@(x,y) (2.*y-1).*(x.^2-x);
fuy_x=@(x,y) (2.*x-1).*(y.^2-y);
fuy_y=@(x,y) (2.*y-1).*(x.^2-x);

% quadrature rule
n_int_xi  = 3;
n_int_eta = 3;
n_int     = n_int_xi * n_int_eta;
[xi, eta, weight] = Gauss2D(n_int_xi, n_int_eta);

% mesh generation
n_en   = 4;               % number of nodes in an element
%different hh
n_el_x0=[20,30,40,50,60,70,80,90]; %changable
n_el_y0=[20,30,40,50,60,70,80,90]; %changable
flag=0;
for n_el_x = n_el_x0              % number of elements in x-dir
    for n_el_y = n_el_y0               % number of elements in y-dir
        n_el   = n_el_x * n_el_y; % total number of elements

        n_np_x = n_el_x + 1;      % number of nodal points in x-dir
        n_np_y = n_el_y + 1;      % number of nodal points in y-dir
        n_np   = n_np_x * n_np_y; % total number of nodal points

        x_coor = zeros(n_np, 1);
        y_coor = x_coor;

        hx = 1.0 / n_el_x;        % mesh size in x-dir
        hy = 1.0 / n_el_y;        % mesh size in y-dir
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
        %have 'flag' different hhs

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
        for ex = 1 : n_el_x
            for ey = 1 : n_el_y
                ee = (ey-1) * n_el_x + ex; % element index
                IEN(ee, 1) = (ey-1) * n_np_x + ex;
                IEN(ee, 2) = (ey-1) * n_np_x + ex + 1;
                IEN(ee, 3) =  ey    * n_np_x + ex + 1;
                IEN(ee, 4) =  ey    * n_np_x + ex;
            end
        end

        boundary_conditions=[0,0,0,0,0,0,0,0];

        disp_top_f = @(x,y,dof) (dof==1)* 0 + (dof==2)* 0;
        disp_bottom_f = @(x,y,dof) (dof==1)* 0+ (dof==2)* 0;
        disp_left_f = @(x,y,dof) (dof==1)* 0 + (dof==2)* 0;
        disp_right_f = @(x,y,dof) (dof==1)* 0 + (dof==2)* 0;

        % ID array
        n_ed=2;
        ID = zeros(n_np,n_ed);
        counter = 0;
        for ny = 2 : n_np_y - 1
            for nx = 2 : n_np_x - 1
                index = (ny-1)*n_np_x + nx;
                counter = counter + 1;
                ID(index,1) = counter;
                counter=counter+1;
                ID(index,2) = counter;
            end
        end

        n_eq = counter;
        LM=zeros(2,n_el,n_en);
        for i=1:n_el
            for j=1:n_en
                LM(1,i,j)=ID( IEN(i,j),1 );
                LM(2,i,j)=ID( IEN(i,j),2 );
            end
        end

        % allocate the stiffness matrix and load vector
        K = spalloc(n_eq, n_eq, 9 * n_eq);
        F = zeros(n_eq, 1);

        % loop over element to assembly the matrix and vector
        for ee = 1 : n_el

            x_ele = x_coor( IEN(ee, 1:n_en) );
            y_ele = y_coor( IEN(ee, 1:n_en) );

            k_ele = zeros(n_ed*n_en,n_ed*n_en); %8*8
            f_ele = zeros(n_ed*n_en,1); %8*1

            for ll = 1 : n_int
                x_l = 0.0; y_l = 0.0;
                dx_dxi = 0.0; dx_deta = 0.0;
                dy_dxi = 0.0; dy_deta = 0.0;

                for aa = 1 : n_en
                    x_l = x_l + x_ele(aa) * Quad(aa, xi(ll), eta(ll));
                    y_l = y_l + y_ele(aa) * Quad(aa, xi(ll), eta(ll));
                    [Na_xi, Na_eta] = Quad_grad(aa, xi(ll), eta(ll));
                    dx_dxi  = dx_dxi  + x_ele(aa) * Na_xi;
                    dx_deta = dx_deta + x_ele(aa) * Na_eta;
                    dy_dxi  = dy_dxi  + y_ele(aa) * Na_xi;
                    dy_deta = dy_deta + y_ele(aa) * Na_eta;
                end

                detJ = dx_dxi * dy_deta - dx_deta * dy_dxi;

                %implement 2
                for aa = 1 : n_en
                    Na = Quad(aa, xi(ll), eta(ll));
                    [Na_xi, Na_eta] = Quad_grad(aa, xi(ll), eta(ll));
                    Na_x = (Na_xi * dy_deta - Na_eta * dy_dxi) / detJ;
                    Na_y = (-Na_xi * dx_deta + Na_eta * dx_dxi) / detJ;

                    Ba=[Na_x, 0; 0, Na_y; Na_y, Na_x];
                    DB=D0*Ba;
                    % B1a=Na_x;B2a=Na_y;
                    % DB(1,1)=D(1,1)*B1a; DB(1,2)=D(1,2)*B2a; DB(2,1)=D(1,2)*B1a; DB(2,2)=D(2,2)*B2a; DB(3,1)=D(3,3)*B2a; DB(3,2)=D(3,3)*B1a;

                    f_ele(2*aa-1) = f_ele(2*aa-1) + weight(ll) * detJ * f(x_l, y_l,1) * Na ;
                    f_ele(2*aa) = f_ele(2*aa) + weight(ll) * detJ * f(x_l, y_l,2) * Na ;
                    for bb = 1 : n_en
                        Nb = Quad(bb, xi(ll), eta(ll));
                        [Nb_xi, Nb_eta] = Quad_grad(bb, xi(ll), eta(ll));
                        Nb_x = (Nb_xi * dy_deta - Nb_eta * dy_dxi) / detJ;
                        Nb_y = (-Nb_xi * dx_deta + Nb_eta * dx_dxi) / detJ;
                        Bb=[Nb_x, 0; 0, Nb_y; Nb_y, Nb_x];
                        BDB=(Bb')*DB;
                        % B1b=Nb_x; B2b=Nb_y;
                        % BDB(1,1)=B1b*DB(1,1)+B2b*DB(3,1); BDB(1,2)=B1b*DB(1,2)+B2b*DB(3,2); BDB(2,1)=B2b*DB(2,1)+B1b*DB(3,1); BDB(2,2)=B2b*DB(2,2)+B1b*DB(3,2);
                        k_ele(2*aa-1,2*bb-1) = k_ele(2*aa-1,2*bb-1) + weight(ll) * detJ * BDB(1,1) ;
                        k_ele(2*aa-1,2*bb) = k_ele(2*aa-1,2*bb)+ weight(ll) * detJ * BDB(1,2) ;
                        k_ele(2*aa,2*bb-1) = k_ele(2*aa,2*bb-1)+ weight(ll) * detJ * BDB(2,1) ;
                        k_ele(2*aa,2*bb) = k_ele(2*aa,2*bb)+ weight(ll) * detJ * BDB(2,2) ;
                    end
                end
            end

            %K and F
            for i=1:2
                for aa = 1 : n_en
                    PP = LM(i,ee, aa);
                    if PP > 0
                        F(PP) = F(PP) + f_ele(2*(aa-1)+i); %p=2*(aa-1)+i q=2*(bb-1)+j
                        for j=1:2
                            for bb = 1 : n_en
                                QQ = LM(j,ee, bb);
                                if QQ > 0
                                    K(PP, QQ) = K(PP, QQ) + k_ele(2*(aa-1)+i, 2*(bb-1)+j);
                                else %QQ=0
                                    %modify F with the boundary data(Dirichlet条件)
                                    %此时都为0，可以不用考虑
                                end
                            end
                        end
                    end
                end
            end
        end

        % solve the stiffness matrix
        dn = K \ F;

        % insert dn back into the vector for all nodes
        disp = zeros(n_np, 2);

        for ii = 1 : n_np
            for i=1:2
                index = ID(ii,i);
                if index > 0
                    disp(ii,i) = dn(index);
                else
                    % modify disp with the g data.此时都为0
                end
            end
        end

        fenzi0_x=0; fenmu0_x=0;
        fenzi1_x=0; fenmu1_x=0;
        fenzi0_y=0; fenmu0_y=0;
        fenzi1_y=0; fenmu1_y=0;
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
                fenmu0_x=
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
    end
end
%plot
%x方向
figure()
plot(log(hh),log(error0_x),'ro');
hold on
A_x(1,:)=polyfit(log(hh),log(error0_x),1);
xa=linspace(-5,-2,20);ya=A_x(1,1)*xa+A_x(1,2);
plot(xa,ya,'-');

plot(log(hh),log(error1_x),'x');
hold on
B_x(1,:)=polyfit(log(hh),log(error1_x),1);
xb=linspace(-5,-2,20);yb=B_x(1,1)*xb+B_x(1,2);
plot(xb,yb,'-')

title('error for x direction');
xlabel('log(h)');
ylabel('log(error)');

%y方向
figure()
plot(log(hh),log(error0_y),'ro');
hold on
A_y(1,:)=polyfit(log(hh),log(error0_y),1);
xa=linspace(-5,-2,20);ya=A_y(1,1)*xa+A_y(1,2);
plot(xa,ya,'-');

plot(log(hh),log(error1_y),'x');
hold on
B_y(1,:)=polyfit(log(hh),log(error1_y),1);
xb=linspace(-5,-2,20);yb=B_y(1,1)*xb+B_y(1,2);
plot(xb,yb,'-')

title('error for y direction');
xlabel('log(h)');
ylabel('log(error)');