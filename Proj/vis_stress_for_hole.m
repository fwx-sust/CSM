%有孔方形的应力可视化
x_plot=zeros(n_el,1);y_plot=zeros(n_el,1);
stress_xx_plot=zeros(n_el,1);stress_yy_plot=zeros(n_el,1);stress_xy_plot=zeros(n_el,1);
stress=zeros(3,1); stress_zz_plot=zeros(n_el,1);
for ee=1:n_el
    x_ele = x_coor( IEN(ee, 1:n_en) );
    y_ele = y_coor( IEN(ee, 1:n_en) );
    ux_ele = disp( IEN(ee, :),1 );
    uy_ele = disp( IEN(ee, :),2 );

    stress_zz=0;

    for qua=1:n_int
        x=0;y=0;
        ux=0;
        uy=0;
        strain_xx=0;
        strain_yy=0;
        strain_xy=0;

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

            x=x+x_ele(aa)*Quad(aa, xi(qua), eta(qua));
            y=y+y_ele(aa)*Quad(aa,xi(qua),eta(qua));
            ux=ux+ux_ele*Quad(aa,xi(qua),eta(qua));
            uy=uy+uy_ele*Quad(aa,xi(qua),eta(qua));

            [Na_xi, Na_eta] = Quad_grad(aa, xi(qua), eta(qua));
            strain_xx=strain_xx+ux_ele(aa)*(Na_xi * dy_deta - Na_eta * dy_dxi) / detJ;
            strain_yy=strain_yy+uy_ele(aa)*(-Na_xi * dx_deta + Na_eta * dx_dxi) / detJ ;
            strain_xy=strain_xy+0.5*ux_ele(aa)*(-Na_xi * dx_deta + Na_eta * dx_dxi) / detJ + 0.5*+uy_ele(aa)*(Na_xi * dy_deta - Na_eta * dy_dxi) / detJ;
        end
    end

    stress=D0*[ strain_xx;strain_yy;2*strain_xy ];

    if question_def == 1
        stress_zz = v*(stress(1,1)+stress(2,1));
        stress_zz_plot(ee)=stress_zz;
    end

    x_plot(ee)=x;y_plot(ee)=y;
    stress_xx_plot(ee)=stress(1,1);stress_yy_plot(ee)=stress(2,1);stress_xy_plot(ee)=stress(3,1);
end

%xx
figure;
% 生成规则网格
[Xq, Yq] = meshgrid(linspace(min(x_plot), max(x_plot), n_el), linspace(min(x_plot), max(y_plot), n_el));

% 插值z到规则网格上
Zq = griddata(x_plot, y_plot,stress_xx_plot, Xq, Yq, 'v4');

% 绘制云图
surf(Xq,Yq,Zq)
shading interp
colorbar;
%view([90, 90]);     % 调整试图位置
hold on;
title('Contour Plot of Function stress for xxh');
xlabel('x');
ylabel('y');

%yy
figure;
% 生成规则网格
[Xq, Yq] = meshgrid(linspace(min(x_plot), max(x_plot), n_el), linspace(min(x_plot), max(y_plot), n_el));

% 插值z到规则网格上
Zq = griddata(x_plot, y_plot,stress_yy_plot, Xq, Yq, 'v4');

% 绘制云图
surf(Xq,Yq,Zq)
shading interp
colorbar;
%view([90, 90]);     % 调整试图位置
hold on;
title('Contour Plot of Function stress for yyh');
xlabel('x');
ylabel('y');

%xy
figure;
% 生成规则网格
[Xq, Yq] = meshgrid(linspace(min(x_plot), max(x_plot), n_el), linspace(min(x_plot), max(y_plot), n_el));

% 插值z到规则网格上
Zq = griddata(x_plot, y_plot,stress_xy_plot, Xq, Yq, 'v4');

% 绘制云图
surf(Xq,Yq,Zq)
shading interp
colorbar;
%view([90, 90]);     % 调整试图位置
hold on;
title('Contour Plot of Function stress for xyh');
xlabel('x');
ylabel('y');

%zz
if question_def== 1 %zz
    figure;
    % 生成规则网格
    [Xq, Yq] = meshgrid(linspace(min(x_plot), max(x_plot), n_el), linspace(min(x_plot), max(y_plot), n_el));

    % 插值z到规则网格上
    Zq = griddata(x_plot, y_plot,stress_zz_plot, Xq, Yq, 'v4');

    % 绘制云图
    surf(Xq,Yq,Zq)
    shading interp
    colorbar;
    %view([90, 90]);     % 调整试图位置
    hold on;
    title('Contour Plot of Function stress for zzh');
    xlabel('x');
    ylabel('y');
end



%exact solution
%xx
figure;
x_plot=zeros(n_np,1);y_plot=zeros(n_np,1);
stress_xx=zeros(n_np,1);stress_yy=zeros(n_np,1);stress_xy=zeros(n_np,1);
stress=zeros(3,1); stress_zz=zeros(n_el,1);
for i=1:n_np
    x_plot(i)=x_coor(i);
    y_plot(i)=y_coor(i);
    stress_xx(i)=f_stress_rr( r(i),theta(i) )*( cos_theta(i) )^2 + f_stress_tt( r(i),theta(i) )*( sin_theta(i) )^2 + f_stress_rt( r(i),theta(i) )*2*cos_theta(i)*sin_theta(i);
    stress_yy(i)=f_stress_rr( r(i),theta(i) )*( sin_theta(i) )^2 + f_stress_tt( r(i),theta(i) )*( cos_theta(i) )^2 - f_stress_rt( r(i),theta(i) )*2*cos_theta(i)*sin_theta(i);
    stress_xy(i)=( f_stress_tt( r(i),theta(i) )-f_stress_rr( r(i),theta(i) ) )*cos_theta(i)*sin_theta(i) + f_stress_rt( r(i),theta(i) )*( ( cos_theta(i) )^2-(sin_theta(i))^2 );
    stress_zz(i)=v*(stress_xx(i)+stress_yy(i));
end
% 生成规则网格
[xq, yq] = meshgrid(linspace(min(x_plot), max(x_plot), n_np), linspace(min(x_plot), max(y_plot), n_np));

% 插值z到规则网格上
zq = griddata(x_plot, y_plot, stress_xx, xq, yq, 'cubic');
surf(xq,yq,zq)
shading interp
colorbar;
%view([90, 90]);     % 调整试图位置
hold on;
title('Contour Plot of Function stress for xx');
xlabel('x');
ylabel('y');

%yy
figure;
% 插值z到规则网格上
zq = griddata(x_plot, y_plot, stress_yy, xq, yq, 'cubic');
surf(xq,yq,zq)
shading interp
colorbar;
%view([90, 90]);     % 调整试图位置
hold on;
title('Contour Plot of Function stress for yy');
xlabel('x');
ylabel('y');

%xy
figure;
% 插值z到规则网格上
zq = griddata(x_plot, y_plot, stress_xy, xq, yq, 'cubic');
surf(xq,yq,zq)
shading interp
colorbar;
%view([90, 90]);     % 调整试图位置
hold on;
title('Contour Plot of Function stress for xy');
xlabel('x');
ylabel('y');

%zz
if question_def ==1 %%zz
    figure;
    zq = griddata(x_plot, y_plot, stress_zz, xq, yq, 'cubic');
    surf(xq,yq,zq)
    shading interp
    colorbar;
    %view([90, 90]);     % 调整试图位置
    hold on;
    title('Contour Plot of Function strain for zz');
    xlabel('x');
    ylabel('y');
end