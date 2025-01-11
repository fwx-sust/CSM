%1*1
%只用来计算四边固支的情况，不可改变，只用来验证内部计算代码
clear all; clc;

E=2e11; %Young's modulus
v=0.3; %Poisson's ratio

%D(x)————specify either plane strain or plane stress
question_def = 1; %1 for plane strain, 2 for plane stress
D0 = stiffnessD(question_def, E, v); %stress=D*strain

%body force
%只适用于plane strain，因为应变和应力的关系为plane strain提供的
Da=E/(v+1)/(2*v-1);
f=@(x,y,dof) (dof==1)*(-Da)*( (v-1)*2*(y^2-y) +(2*v-1)*(x^2-x)*2 ) +...
    (dof==2)*(-Da)*( (2*v-1)*(2*x-1)*(2*y-1) );

%exact solution————（假设为二次位移场） 以下是最简单的两端固支的位移

fux=@(x,y) x.*(x-1).*y.*(y-1); %点乘方便后续画图
fuy=@(x,y) 0;

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

%求单元面积
area=zeros(n_el);
area_s=zeros(n_el);
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

boundary_conditions=[0,0,0,0,0,0,0,0];

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
        D=weight(ll)*detJ*D0; %set up D~

        for aa = 1 : n_en
            Na = Quad(aa, xi(ll), eta(ll));
            [Na_xi, Na_eta] = Quad_grad(aa, xi(ll), eta(ll));
            Na_x = (Na_xi * dy_deta - Na_eta * dy_dxi) / detJ;
            Na_y = (-Na_xi * dx_deta + Na_eta * dx_dxi) / detJ;
            
            Ba=[Na_x, 0; 0, Na_y; Na_y, Na_x];
            DB=D*Ba;
            % B1a=Na_x;B2a=Na_y;
            % DB(1,1)=D(1,1)*B1a; DB(1,2)=D(1,2)*B2a; DB(2,1)=D(1,2)*B1a; DB(2,2)=D(2,2)*B2a; DB(3,1)=D(3,3)*B2a; DB(3,2)=D(3,3)*B1a;

            f_ele(2*aa-1) = f_ele(2*aa-1) + weight(ll) * detJ * f(x_l, y_l,1) * Na;
            f_ele(2*aa) = f_ele(2*aa) + weight(ll) * detJ * f(x_l, y_l,2) * Na;
            for bb = 1 : aa
                Nb = Quad(bb, xi(ll), eta(ll));
                [Nb_xi, Nb_eta] = Quad_grad(bb, xi(ll), eta(ll));
                Nb_x = (Nb_xi * dy_deta - Nb_eta * dy_dxi) / detJ;
                Nb_y = (-Nb_xi * dx_deta + Nb_eta * dx_dxi) / detJ;
                Bb=[Nb_x, 0; 0, Nb_y; Nb_y, Nb_x];
                BDB=(Bb')*DB;
                % B1b=Nb_x; B2b=Nb_y;
                % BDB(1,1)=B1b*DB(1,1)+B2b*DB(3,1); BDB(1,2)=B1b*DB(1,2)+B2b*DB(3,2); BDB(2,1)=B2b*DB(2,1)+B1b*DB(3,1); BDB(2,2)=B2b*DB(2,2)+B1b*DB(3,2);
                k_ele(2*aa-1,2*bb-1) = k_ele(2*aa-1,2*bb-1) + weight(ll) * detJ * BDB(1,1) * (Na_x * Nb_x + Na_y * Nb_y);
                k_ele(2*aa-1,2*bb) = k_ele(2*aa-1,2*bb)+ weight(ll) * detJ * BDB(1,2) * (Na_x * Nb_x + Na_y * Nb_y);
                k_ele(2*aa,2*bb-1) = k_ele(2*aa,2*bb-1)+ weight(ll) * detJ * BDB(2,1) * (Na_x * Nb_x + Na_y * Nb_y);
                k_ele(2*aa,2*bb) = k_ele(2*aa,2*bb)+ weight(ll) * detJ * BDB(2,2) * (Na_x * Nb_x + Na_y * Nb_y);
            end
        end
    end

    %K and F
    for aa = 1 : n_en
        for i=1:2
            PP = LM(i,ee, aa);
            if PP > 0
                F(PP) = F(PP) + f_ele(2*(aa-1)+i); %p=2*(aa-1)+i q=2*(bb-1)+i
                for bb = 1 : n_en
                    QQ = LM(i,ee, bb);
                    if QQ > 0
                        K(PP, QQ) = K(PP, QQ) + k_ele(2*(aa-1)+i, 2*(bb-1)+i);
                    else %QQ=0
                        % modify F with the boundary data(Dirichlet条件)
                        % 此时都为0，可以不用考虑
                        % a=x_coor( IEN(ee,bb) );
                        % b=y_coor( IEN(ee,bb) );
                        % for j=1:size(top_pos,1)
                        %     if a== top_pos(j,1) && b==top_pos(j,2)
                        %         F(PP)= F(PP) - k_ele(2*(aa-1)+i, 2*(bb-1)+i)*disp_top_f(a,b,i);
                        %     end
                        % end
                        % for j=1:size(bottom_pos,1)
                        %     if a== bottom_pos(j,1) && b==bottom_pos(j,2)
                        %         F(PP)= F(PP) - k_ele(2*(aa-1)+i, 2*(bb-1)+i)*disp_bottom_f(a,b,i);
                        %     end
                        % end
                        % for j=1:size(left_pos,1)
                        %     if a== left_pos(j,1) && b==left_pos(j,2)
                        %         F(PP)= F(PP) - k_ele(2*(aa-1)+i, 2*(bb-1)+i)*disp_left_f(a,b,i);
                        %     end
                        % end
                        % for j=1:size(right_pos,1)
                        %     if a== right_pos(j,1) && b==right_pos(j,2)
                        %         F(PP)= F(PP) - k_ele(2*(aa-1)+i, 2*(bb-1)+i)*disp_right_f(a,b,i);
                        %     end
                        % end
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
            a=x_coor( ii );
            b=y_coor( ii );
            for j=1:size(top_pos,1)
                if a== top_pos(j,1) && b==top_pos(j,2)
                    disp(ii,i) = disp_top_f(a,b,i);
                end
            end
            for j=1:size(bottom_pos,1)
                if a== bottom_pos(j,1) && b==bottom_pos(j,2)
                    disp(ii,i) = disp_bottom_f(a,b,i);
                end
            end
            for j=1:size(left_pos,1)
                if a== left_pos(j,1) && b==left_pos(j,2)
                    disp(ii,i) = disp_left_f(a,b,i);
                end
            end
            for j=1:size(right_pos,1)
                if a== right_pos(j,1) && b==right_pos(j,2)
                    disp(ii,i) = disp_right_f(a,b,i);
                end
            end
        end
    end
end
%以上就算出来了位移disp

%整理得到uhx,uhy，strain,stress
x_plot=[]; y_plot=[];uhx_plot=[];uhy_plot=[];
for ee=1:n_el
    x_ele = x_coor( IEN(ee, :) );
    y_ele = y_coor( IEN(ee ,:) );
    u_ele_x = disp( IEN(ee, :),1 );
    u_ele_y = disp( IEN(ee, :),2 );
    epsilon=zeros(3,1);
    stress=zeros(3,1);
    for qua=1:n_int
        uhx=0;
        uhy=0;
        uhx_x=0;uhx_y=0;
        uhy_x=0;uhy_y=0;
        x=0;
        y=0;
        for aa = 1 : n_en
            x=x+x_ele(aa)*Quad(aa, xi(qua), eta(qua));
            y=y+y_ele(aa)*Quad(aa,xi(qua),eta(qua));

            %disp
            uhx=uhx+u_ele_x(aa)*Quad(aa, xi(qua), eta(qua));
            uhy=uhy+u_ele_y(aa)*Quad(aa, xi(qua), eta(qua));

            [Na_xi, Na_eta] = Quad_grad(aa, xi(qua), eta(qua));

            uhx_x=uhx_x+u_ele_x(aa)*(Na_xi * dy_deta - Na_eta * dy_dxi) / detJ;
            uhx_y=uhx_y+u_ele_x(aa)*(-Na_xi * dx_deta + Na_eta * dx_dxi) / detJ;

            uhy_x=uhx_y+u_ele_y(aa)*(Na_xi * dy_deta - Na_eta * dy_dxi) / detJ;
            uhy_y=uhy_y+u_ele_y(aa)*(-Na_xi * dx_deta + Na_eta * dx_dxi) / detJ ;

            %strain
            epsilon(1,1)=epsilon(1,1)+uhx_x; %11
            epsilon(2,1)=epsilon(2,1)+uhy_y; %22
            epsilon(3,1)=epsilon(3,1)+0.5*uhx_y+0.5*uhy_x;%12

            %stress
            stress=D*epsilon;
        end
    end
    x_plot=[x_plot,x];
    y_plot=[y_plot,y];
    uhx_plot=[uhx_plot,uhx];
    uhy_plot=[uhy_plot,uhy];
end

%绘制displacement
% 创建一个规则网格
[xq, yq] = meshgrid(linspace(min(x_plot), max(x_plot), 100),linspace(min(y_plot), max(y_plot), 100));

% 使用 griddata 函数进行插值
zq = griddata(x_plot, y_plot, uhy_plot, xq, yq, 'cubic');

% 绘制云图
figure;
contourf(xq, yq, zq, 20, 'LineColor', 'none');
colorbar; % 显示颜色条
hold on;
title('Contour Plot of Function uhy');
xlabel('x');
ylabel('y');

%画的y
figure;
plot3(xq,yq,zq)
xlabel("x轴")
ylabel("y轴")
zlabel("z轴")
grid on

figure;
zq= griddata(x_plot, y_plot, uhx_plot, xq, yq, 'cubic');
contourf(xq, yq, zq, 20, 'LineColor', 'none');
colorbar; % 显示颜色条
hold on;
title('Contour Plot of Function uhx');
xlabel('x');
ylabel('y')


%画的x
figure;
plot3(xq,yq,zq)
xlabel("x轴")
ylabel("y轴")
zlabel("z轴")
grid on

%exact solution————上边界受拉力，下边界固定
figure;

x = linspace(min(x_plot), max(x_plot), 50);
y = linspace(min(x_plot), max(x_plot), 50);
[xq,yq]=meshgrid(x,y);
uy=zeros(50,50); %和x,y所分份数有关
uy(:,:)=fuy(xq,yq);

zq = griddata(x, y, uy, xq, yq, 'cubic'); % 使用三次插值
contourf(xq, yq, zq, 20, 'LineColor', 'none');
colorbar; % 显示颜色条
hold on;
title('Contour Plot of Function uy');
xlabel('x');
ylabel('y');

figure;
ux=zeros(50,50);
ux(:,:)=fux(xq,yq);
zq = griddata(x, y, ux, xq, yq, 'cubic'); % 使用三次插值
contourf(xq, yq, zq, 20, 'LineColor', 'none');
colorbar; % 显示颜色条
hold on;
title('Contour Plot of Function ux');
xlabel('x');
ylabel('y');

%绘制三维图像
% 给定x、y、z的数值
figure;
plot3(xq,yq,zq)
xlabel("x轴")
ylabel("y轴")
zlabel("z轴")
grid on