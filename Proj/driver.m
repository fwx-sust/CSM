%drive
%1*1矩形实验
clear all; clc;

E=2e11; %Young's modulus
v=0.3; %Poisson's ratio
q=1e5; %N/m^2

%D(x)————specify either plane strain or plane stress
question_def = 2; %1 for plane strain, 2 for plane stress
D = stiffnessD(question_def, E, v);

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
end

%力边界条件: 可改
f=@(x,y,dof) 0;%no body force

% =0 代表Dirichlet（设置位移g） Neumann: =1 代表应力边界条件（仅考虑自由边界，及应力为0的情况） =2代表设置力边界条件
%可变
boundary_conditions=[2,2,1,1,1,1,0,0]; %顺序为上下左右(各有两个方向） 下边界固定，左右边界自由，上边界受力

%提供Dirichlet条件g————位移边界条件：可以改变
% for i_bc=1:8 %四条边*2
%     if boundary_conditions(i_bc) == 0
%         if i_bc==1 || i_bc==2
%             %(dof == 1: x) (dof == 2: y)
%             disp_top_f = @(x,y,dof) (dof==1)* 0 + (dof==2)* 0;
%         else
%             if  i_bc==3 || i_bc==4
%                 disp_bottom_f = @(x,y,dof) (dof==1)* 0 + (dof==2)* 0;
%             else
%                 if i_bc==5 || i_bc==6
%                     disp_left_f = @(x,y,dof) (dof==1)* 0 + (dof==2)* 0;
%                 else
%                     disp_right_f = @(x,y,dof) (dof==1)* 0 + (dof==2)* 0;
%                 end
%             end
%         end
%     end
% end

%手动输入
disp_top_f = @(x,y,dof) (dof==1)* 0 + (dof==2)* 0;
%除了下边界，其余均是Not necessary for the question I defined
disp_bottom_f = @(x,y,dof) (dof==1)* 0 + (dof==2)* 0;
disp_left_f = @(x,y,dof) (dof==1)* 0 + (dof==2)* 0;
disp_right_f = @(x,y,dof) (dof==1)* 0 + (dof==2)* 0;

%提供Neumann边界条件h————应力边界条件：只考虑应力为0的情况，不可改
% for i_bc=1:8
%     if boundary_conditions(i_bc) == 1
%         if i_bc==1
%             %(dof == 1: x) (dof == 2: y)
%             stress_top_f = @(x,y,dof) (dof==1)* 0 + (dof==2)* 0;
%         else
%             if  i_bc==2
%                 stress_bottom_f = @(x,y,dof) (dof==1)* 0 + (dof==2)* 0;
%             else
%                 if i_bc==3
%                     stress_left_f = @(x,y,dof) (dof==1)* 0 + (dof==2)* 0;
%                 else
%                     stress_right_f = @(x,y,dof) (dof==1)* 0 + (dof==2)* 0;
%                 end
%             end
%         end
%     end
% end

%手动输入，但该代码只考虑自由端，即应力为0
%同理，除了左右边界，其余的not necessary
stress_top_f = @(x,y,dof) (dof==1)* 0 + (dof==2)* 0;
stress_bottom_f = @(x,y,dof) (dof==1)* 0 + (dof==2)* 0;
stress_left_f = @(x,y,dof) (dof==1)* 0 + (dof==2)* 0;
stress_right_f = @(x,y,dof) (dof==1)* 0 + (dof==2)* 0;

% ID array——两个自由度
ID = zeros(n_np,n_ed);
counter = 0;
A_for_q=[];
for i=1:n_np
    %判断是否在边界上
    flag=0; %=1 top; =2 bottom; =3 left; =4 right 和bc矩阵对应
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
        A_for_q(counter) =0;
        counter=counter+1;
        ID(i,2)=counter;
        A_for_q(counter) =0;
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
        %记录外载荷向量的位置
        if boundary_conditions(flag+1) ==2
            A_for_q(counter)=1; %y方向
            A_for_q(counter-1) =0;
        else
            A_for_q(counter)=0;
            A_for_q(counter-1)=0;
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

        for aa = 1 : n_en
            Na = Quad(aa, xi(ll), eta(ll));
            [Na_xi, Na_eta] = Quad_grad(aa, xi(ll), eta(ll));
            Na_x = (Na_xi * dy_deta - Na_eta * dy_dxi) / detJ;
            Na_y = (-Na_xi * dx_deta + Na_eta * dx_dxi) / detJ;

            Ba=[Na_x, 0; 0, Na_y; Na_y, Na_x];
            kappa1=D*Ba;

            f_ele(2*aa-1) = f_ele(2*aa-1) + weight(ll) * detJ * f(x_l, y_l,1) * Na;
            f_ele(2*aa) = f_ele(2*aa) + weight(ll) * detJ * f(x_l, y_l,2) * Na;
            for bb = 1 : aa
                Nb = Quad(bb, xi(ll), eta(ll));
                [Nb_xi, Nb_eta] = Quad_grad(bb, xi(ll), eta(ll));
                Nb_x = (Nb_xi * dy_deta - Nb_eta * dy_dxi) / detJ;
                Nb_y = (-Nb_xi * dx_deta + Nb_eta * dx_dxi) / detJ;
                Bb=[Nb_x, 0; 0, Nb_y; Nb_y, Nb_x];
                kappa=(Bb')*kappa1;
                k_ele(2*aa-1,2*bb-1) = k_ele(2*aa-1,2*bb-1) + weight(ll) * detJ * kappa(1,1) * (Na_x * Nb_x + Na_y * Nb_y);
                k_ele(2*aa-1,2*bb) = k_ele(2*aa-1,2*bb)+ weight(ll) * detJ * kappa(1,2) * (Na_x * Nb_x + Na_y * Nb_y);
                k_ele(2*aa,2*bb-1) = k_ele(2*aa,2*bb-1)+ weight(ll) * detJ * kappa(2,1) * (Na_x * Nb_x + Na_y * Nb_y);
                k_ele(2*aa,2*bb) = k_ele(2*aa,2*bb)+ weight(ll) * detJ * kappa(2,2) * (Na_x * Nb_x + Na_y * Nb_y);
            end
        end
    end

    %K and F

    for aa = 1 : n_en
        for i=1:2
            PP = LM(i,ee, aa);
            if PP > 0
                F(PP) = F(PP) + f_ele(2*(aa-1)+i);
                for bb = 1 : n_en
                    QQ = LM(i,ee, bb);
                    if QQ > 0
                        K(PP, QQ) = K(PP, QQ) + k_ele(2*(aa-1)+i, 2*(bb-1)+i);
                    else %QQ=0
                        % modify F with the boundary data(Dirichlet条件)
                        a=x_coor( IEN(ee,bb) );
                        b=y_coor( IEN(ee,bb) );
                        for j=1:size(top_pos,1)
                            if a== top_pos(j,1) && b==top_pos(j,2)
                                F(PP)= F(PP) - k_ele(2*(aa-1)+i, 2*(bb-1)+i)*disp_top_f(a,b,i);
                            end
                        end
                        for j=1:size(bottom_pos,1)
                            if a== bottom_pos(j,1) && b==bottom_pos(j,2)
                                F(PP)= F(PP) - k_ele(2*(aa-1)+i, 2*(bb-1)+i)*disp_bottom_f(a,b,i);
                            end
                        end
                        for j=1:size(left_pos,1)
                            if a== left_pos(j,1) && b==left_pos(j,2)
                                F(PP)= F(PP) - k_ele(2*(aa-1)+i, 2*(bb-1)+i)*disp_left_f(a,b,i);
                            end
                        end
                        for j=1:size(right_pos,1)
                            if a== right_pos(j,1) && b==right_pos(j,2)
                                F(PP)= F(PP) - k_ele(2*(aa-1)+i, 2*(bb-1)+i)*disp_right_f(a,b,i);
                            end
                        end
                    end
                end
            end
        end
    end
    %在此处考虑添加的力边界条件，修改F
    for i_bc=1:4
        if boundary_conditions(i_bc) == 2
            %但是我不会
        end
    end
end
%如果只在端点处加力边界条件，在此处修改F
F=q*A_for_q';

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