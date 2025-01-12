%drive
%square with hole
%Fig.1 
%极坐标系
clear all; clc;

%自定义
E=1e9; %Young's modulus
v=0.3; %Poisson's ratio

%body force: 从应力分布推导
f=@(x,y,dof) (dof==1)*0 + (dof==2)*0; 

%exact solution
%应力理论值
Tx=10*1e3; %N/m^2 Pa
R=0.5;
f_stress_rr = @(r,theta) 0.5*Tx.*( 1-(R^2)./(r.^2))+0.5*Tx*(1-4*(R^2)./(r.^2)+3*(R^4)./(r.^4))*cos(2.*theta);
f_stress_tt = @(r,theta) 0.5*Tx.*( 1+(R^2)./(r.^2))-0.5*Tx*(1+3*(R^4)./(r.^4))*cos(2.*theta);
f_stress_rt = @(r,theta) -0.5*Tx*(1+2*(R^2)./(r.^2)-3*(R^4)./(r.^4))*sin(2.*theta);

%D(x)————specify either plane strain or plane stress：自定义
question_def = 1; %1 for plane strain, 2 for plane stress
D0 = stiffnessD(question_def, E, v);

% quadrature rule
n_int_xi  = 3;
n_int_eta = 3;
n_int     = n_int_xi * n_int_eta;
[xi, eta, weight] = Gauss2D(n_int_xi, n_int_eta);

% mesh generation——gmsh生成
%n_np:total number of nodal points
%n_el:total number of elements
%改变网格可能需要改变读网格的方法，但若只是更改每条边的节点数，可以直接用这个方法
[n_np,n_el,x_coor,y_coor,IEN,top_pos,bottom_pos,left_pos,right_pos,circle_pos] = read_mesh_with_hole('quarter_plate_with_hole.m');
%top_pos等返回的是边界坐标，后面发现其实直接返回A会更好

%r and theta
r=zeros(n_np,1);theta=zeros(n_np,1);
sin_theta=zeros(n_np,1); cos_theta=zeros(n_np,1);
for i=1:n_np
    r(i) = sqrt(x_coor(i)^2 + y_coor(i)^2);
    theta(i) = atan( y_coor(i)/x_coor(i) );
    sin_theta(i) = y_coor(i) / r(i);
    cos_theta(i) = x_coor(i) / r(i);
end

n_en = 4;  % number of nodes in an element 四边形网格
n_ed = 2;  %两个自由度

%求单元面积
area=zeros(n_el,1);
area_s=zeros(n_el,1);
for i=1:n_el
    x1=x_coor( IEN(i,1) );
    y1=y_coor( IEN(i,1) );
    x2=x_coor( IEN(i,2) );
    y2=y_coor( IEN(i,2) );
    x3=x_coor( IEN(i,3) );
    y3=y_coor( IEN(i,3) );
    x4=x_coor( IEN(i,4) );
    y4=y_coor( IEN(i,4) );
    area(i)=0.5 * abs( x1*y2 + x2*y3 + x3*y4 + x4*y1 - (y1*x2 + y2*x3 + y3*x4 + y4*x1) );
    area_s(i)=sqrt(area(i));
end
h=max(area_s);

%边界条件

%外载荷边界条件
force_boundary=[0,0,0,0,0]; %==0代表没有均布载荷，==1代表有均布载荷（只考虑均布载荷）
force_top_f = @(x,y,dof) (dof==1)* 0 + (dof==2)* 0;
force_bottom_f = @(x,y,dof) (dof==1)* 0+ (dof==2)* 0;
force_left_f = @(x,y,dof) (dof==1)* 0 + (dof==2)* 0;
force_right_f = @(x,y,dof) (dof==1)* Tx + (dof==2)* 0;
force_circle_f=@(x,y,dof) (dof==1)* 0 + (dof==2)* 0;

% =0 代表Dirichlet（设置位移g） Neumann: =1 代表应力边界条件（设置应力）
%自定义
boundary_conditions=[1,1, 1,0, 0,1, 1,1, 1,1]; %有五条边,上下左右圆
%提供Dirichlet条件g————位移边界条件：可以改变
%手动输入————dof==1是x方向，dof==2是y方向
disp_top_f = @(x,y,dof) (dof==1)* 0 + (dof==2)* 0;
disp_bottom_f = @(x,y,dof) (dof==1)* 0+ (dof==2)* 0;
disp_left_f = @(x,y,dof) (dof==1)* 0 + (dof==2)* 0;
disp_right_f = @(x,y,dof) (dof==1)* 0 + (dof==2)* 0;
disp_circle_f=@(x,y,dof) (dof==1)* 0 + (dof==2)* 0;

%提供Neumann边界条件h————应力边界条件：可以改变
%手动输入————dof==1是x方向，dof==2是y方向
%上，右要改
stress_top_f = @(x,y,dof) (dof==1)* 0 + (dof==2)* 0;
stress_bottom_f = @(x,y,dof) (dof==1)* 0 + (dof==2)* 0;
stress_left_f = @(x,y,dof) (dof==1)* 0 + (dof==2)* 0;
stress_right_f = @(x,y,dof) (dof==1)* 0 + (dof==2)* 0;
stress_circle_f=@(x,y,dof)  (dof==1)* 0 + (dof==2)* 0;

% ID array——两个自由度
ID = zeros(n_np,n_ed);
counter = 0;
for i=1:n_np
    %判断是否在边界上
    flag=0;  %=1 top; =3 bottom; =5 left; 7 right 和bc矩阵对应
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

    for j=1:size(circle_pos,1)
        if a==circle_pos(j,1) && b==circle_pos(j,2)
            flag=9;
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
        if boundary_conditions(flag) ~=0 %x方向
            counter=counter+1;
            ID(i,1)=counter;
        end
        if boundary_conditions(flag+1) ~=0 %y方向
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
K = spalloc(n_eq, n_eq,n_eq^2 -2);
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

        for aa = 1 : n_en
            Na = Quad(aa, xi(ll), eta(ll));
            [Na_xi, Na_eta] = Quad_grad(aa, xi(ll), eta(ll));
            Na_x = (Na_xi * dy_deta - Na_eta * dy_dxi) / detJ;
            Na_y = (-Na_xi * dx_deta + Na_eta * dx_dxi) / detJ;

            Ba=[Na_x, 0; 0, Na_y; Na_y, Na_x];
            DB=D0*Ba;
            % B1a=Na_x;B2a=Na_y;
            % DB(1,1)=D(1,1)*B1a; DB(1,2)=D(1,2)*B2a; DB(2,1)=D(1,2)*B1a; DB(2,2)=D(2,2)*B2a; DB(3,1)=D(3,3)*B2a; DB(3,2)=D(3,3)*B1a;

            f_ele(2*aa-1) = f_ele(2*aa-1) + weight(ll) * detJ * f(x_l, y_l,1) * Na;
            f_ele(2*aa) = f_ele(2*aa) + weight(ll) * detJ * f(x_l, y_l,2) * Na;

             % modify F with the boundary data(Neumann条件:stress 和均布载荷）
            a=x_coor( IEN(ee,aa) );
            b=y_coor( IEN(ee,aa) );
            %先判断bc，再判断位置能减少计算时间
            if boundary_conditions(1)==1 || boundary_conditions(2)==1
                for m=1:size(top_pos,1)
                    if a== top_pos(m,1) && b==top_pos(m,2)
                        f_ele(2*aa-1) = f_ele(2*aa-1) + weight(ll) * detJ * stress_top_f(a,b,1) * Na;
                        f_ele(2*aa) = f_ele(2*aa) + weight(ll) * detJ * stress_top_f(a,b,2) * Na;
                        if force_boundary(1) == 1   %top一般只受y方向上的均布载荷
                            f_ele(2*aa-1) = f_ele(2*aa-1) + weight(ll)*detJ*force_top_f(a,b,1)*Na;
                            f_ele(2*aa) = f_ele(2*aa) + weight(ll)*detJ*force_top_f(a,b,2)*Na;
                        end
                    end
                end
            end
            if boundary_conditions(3)==1 || boundary_conditions(4)==1
                for m=1:size(bottom_pos,1)
                    if a== bottom_pos(m,1) && b==bottom_pos(m,2)
                        f_ele(2*aa-1) = f_ele(2*aa-1) + weight(ll) * detJ * stress_bottom_f(a,b,1) * Na;
                        f_ele(2*aa) = f_ele(2*aa) + weight(ll) * detJ * stress_bottom_f(a,b,2) * Na;
                        if force_boundary(2) == 1
                            f_ele(2*aa-1) = f_ele(2*aa-1) + weight(ll)*detJ*force_bottom_f(a,b,1)*Na;
                            f_ele(2*aa) = f_ele(2*aa) + weight(ll)*detJ*force_bottom_f(a,b,2)*Na;
                        end
                    end
                end
            end
            if boundary_conditions(5)==1 || boundary_conditions(6)==1
                for m=1:size(left_pos,1)
                    if a== left_pos(m,1) && b==left_pos(m,2)
                        f_ele(2*aa-1) = f_ele(2*aa-1) + weight(ll) * detJ * stress_left_f(a,b,1) * Na;
                        f_ele(2*aa) = f_ele(2*aa) + weight(ll) * detJ * stress_left_f(a,b,2) * Na;
                        if force_boundary(3) == 1
                            f_ele(2*aa-1) = f_ele(2*aa-1) + weight(ll)*detJ*force_left_f(a,b,1)*Na;
                            f_ele(2*aa) = f_ele(2*aa) + weight(ll)*detJ*force_left_f(a,b,2)*Na;
                        end
                    end
                end
            end
            if boundary_conditions(7)==1 || boundary_conditions(8)==1
                for m=1:size(right_pos,1)
                    if a== right_pos(m,1) && b==right_pos(m,2)
                        f_ele(2*aa-1) = f_ele(2*aa-1) + weight(ll) * detJ * stress_right_f(a,b,1) * Na;
                        f_ele(2*aa) = f_ele(2*aa) + weight(ll) * detJ * stress_right_f(a,b,2) * Na;
                        if force_boundary(4) == 1
                            f_ele(2*aa-1) = f_ele(2*aa-1) + weight(ll)*detJ*force_right_f(a,b,1)*Na;
                            f_ele(2*aa) = f_ele(2*aa) + weight(ll)*detJ*force_right_f(a,b,2)*Na;
                        end
                    end
                end
            end

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
    for aa = 1 : n_en
        for i=1:2
            PP = LM(i,ee, aa);
            if PP > 0
                F(PP) = F(PP) + f_ele(2*(aa-1)+i);
                for j=1:2
                    for bb = 1 : n_en
                        QQ = LM(j,ee, bb);
                        if QQ > 0
                            K(PP, QQ) = K(PP, QQ) + k_ele(2*(aa-1)+i, 2*(bb-1)+j);
                        else %QQ=0
                            % modify F with the boundary data(Dirichlet条件)
                            a=x_coor( IEN(ee,bb) );
                            b=y_coor( IEN(ee,bb) );
                            %先判断bc,在判断位置能简化计算
                            if boundary_conditions(j)==0
                                for m=1:size(top_pos,1)
                                    if a== top_pos(m,1) && b==top_pos(m,2)
                                        F(PP)= F(PP) - k_ele(2*(aa-1)+i, 2*(bb-1)+j)*disp_top_f(a,b,j);
                                    end
                                end
                            end
                            if boundary_conditions(2+j)==0
                                for m=1:size(bottom_pos,1)
                                    if a== bottom_pos(m,1) && b==bottom_pos(m,2)
                                        F(PP)= F(PP) - k_ele(2*(aa-1)+i, 2*(bb-1)+j)*disp_bottom_f(a,b,j);
                                    end
                                end
                            end
                            if boundary_conditions(4+j)==0
                                for m=1:size(left_pos,1)
                                    if a== left_pos(m,1) && b==left_pos(m,2)
                                        F(PP)= F(PP) - k_ele(2*(aa-1)+i, 2*(bb-1)+j)*disp_left_f(a,b,j);
                                    end
                                end
                            end
                            if boundary_conditions(6+j)==0
                                for m=1:size(right_pos,1)
                                    if a== right_pos(m,1) && b==right_pos(m,2)
                                        F(PP)= F(PP) - k_ele(2*(aa-1)+i, 2*(bb-1)+j)*disp_right_f(a,b,j);
                                    end
                                end
                            end
                            if boundary_conditions(8+j)==0
                                for m=1:size(circle_pos,1)
                                    if a== circle_pos(m,1) && b==circle_pos(m,2)
                                        F(PP)= F(PP) - k_ele(2*(aa-1)+i, 2*(bb-1)+j)*disp_circle_f(a,b,j);
                                    end
                                end
                            end
                        end
                    end
                end
            else %PP=0

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
            % modify disp with the g data.
            a=x_coor( ii );
            b=y_coor( ii );
            %判断在哪条边上
            if boundary_conditions(i)==0
                for j=1:size(top_pos,1)
                    if a== top_pos(j,1) && b==top_pos(j,2)
                        disp(ii,i) = disp_top_f(a,b,i);
                    end
                end
            end
            if boundary_conditions(2+i)==0
                for j=1:size(bottom_pos,1)
                    if a== bottom_pos(j,1) && b==bottom_pos(j,2)
                        disp(ii,i) = disp_bottom_f(a,b,i);
                    end
                end
            end
            if boundary_conditions(4+i)==0
                for j=1:size(left_pos,1)
                    if a== left_pos(j,1) && b==left_pos(j,2)
                        disp(ii,i) = disp_left_f(a,b,i);
                    end
                end
            end
            if boundary_conditions(6+i)==0
                for j=1:size(right_pos,1)
                    if a== right_pos(j,1) && b==right_pos(j,2)
                        disp(ii,i) = disp_right_f(a,b,i);
                    end
                end
            end
            if boundary_conditions(8+i)==0
                for j=1:size(circle_pos,1)
                    if a== circle_pos(j,1) && b==circle_pos(j,2)
                        disp(ii,i) = disp_circle_f(a,b,i);
                    end
                end
            end
        end
    end
end
%以上就算出来了位移disp