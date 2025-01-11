%1*1
%只用来计算四边固支的情况，不可改变，只用来验证内部计算代码
clear all;

E=100; %Young's modulus
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
%应变理论值
f_strain_xx=@(x,y) (2.*x-1).*(y.^2-y);
f_strain_yy=@(x,y) (x.^2-x).*(2.*y-1);
f_strain_xy=@(x,y) 0.5*( (x.^2-x).*(2.*y-1)+(2.*x-1).*(y.^2-y) );
f_strain_zz=@(x,y) -v/E *( Da*( (v-1).*(2.*x-1).*(y.^2-y)-v.*(x.^2-x).*(2.*y-1) )+ Da*( -v.*(2.*x-1).*(y.^2-y)+(v-1).*(x.^2-x).*(2.*y-1) ) ); %only for plane stress
%应力理论值
f_stress_xx=@(x,y) Da*( (v-1).*(2.*x-1).*(y.^2-y)-v.*(x.^2-x).*(2.*y-1) );
f_stress_yy=@(x,y) Da*( -v.*(2.*x-1).*(y.^2-y)+(v-1).*(x.^2-x).*(2.*y-1) );
f_stress_xy=@(x,y) Da*(2*v-1)*0.5.*( (x.^2-x).*(2.*y-1)+(2.*x-1).*(y.^2-y) );
f_stress_zz=@(x,y) v*( Da*( (v-1).*(2.*x-1).*(y.^2-y)-v.*(x.^2-x).*(2.*y-1) )+Da*( -v.*(2.*x-1).*(y.^2-y)+(v-1).*(x.^2-x).*(2.*y-1) ) ); %only for plane strain

% quadrature rule
n_int_xi  = 3;
n_int_eta = 3;
n_int     = n_int_xi * n_int_eta;
[xi, eta, weight] = Gauss2D(n_int_xi, n_int_eta);

% mesh generation——gmsh生成
%n_np:total number of nodal points
%n_el:total number of elements
[n_np,n_el,x_coor,y_coor,IEN,top_pos,bottom_pos,left_pos,right_pos] = try_to_read_mesh('rectangle_gmsh.m');
n_en = 4;               % number of nodes in an element 四边形网格
n_ed = 2;  %两个自由度

boundary_conditions=[0,0,0,0,0,0,0,0];

%因为全是Dirichlet条件，所以直接设置位移边界条件，不管其他条件
disp_top_f = @(x,y,dof) (dof==1)* 0 + (dof==2)* 0;
disp_bottom_f = @(x,y,dof) (dof==1)* 0+ (dof==2)* 0;
disp_left_f = @(x,y,dof) (dof==1)* 0 + (dof==2)* 0;
disp_right_f = @(x,y,dof) (dof==1)* 0 + (dof==2)* 0;

% ID array——两个自由度
ID = zeros(n_np,n_ed);
counter = 0;
for i=1:n_np
    %判断是否在边界上
    flag=0; %=1 top; =3 bottom; =5 left; 7 right 和bc矩阵对应
    a=x_coor(i);
    b=y_coor(i);

    for j=1:size(top_pos,1)
        if a== top_pos(j,1) && b==top_pos(j,2)
            flag=1;
        end
    end

    for j=1:size(bottom_pos,1)
        if a==bottom_pos(j,1) && b==bottom_pos(j,2)
            flag=3;
        end
    end

    for j=1:size(left_pos,1)
        if a==left_pos(j,1) && b==left_pos(j,2)
            flag=5;
        end
    end

    for j=1:size(right_pos,1)
        if a==right_pos(j,1) && b==right_pos(j,2)
            flag=7;
        end
    end

    %不在边界上则直接赋P值
    if flag==0
        counter=counter+1;
        ID(i,1)=counter;
        counter=counter+1;
        ID(i,2)=counter;
    else
        %在边界上，Neumann边界照常给P，Dirichlet直接默认0
        if boundary_conditions(flag) ~=0
            counter=counter+1;
            ID(i,1)=counter;
        end
        if boundary_conditions(flag+1) ~=0
            counter=counter+1;
            ID(i,2)=counter;
        end
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
K = spalloc(n_eq, n_eq,9*n_eq);
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
%displacement
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
%以上就算出来了位移disp