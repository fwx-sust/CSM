%应变可视化————surf函数
x_plot=zeros(n_el,1);y_plot=zeros(n_el,1);
strain_xx_plot=zeros(n_el,1);strain_yy_plot=zeros(n_el,1);strain_xy_plot=zeros(n_el,1);
stress=zeros(3,1); strain_zz_plot=zeros(n_el,1);
for ee=1:n_el
    x_ele = x_coor( IEN(ee, 1:n_en) );
    y_ele = y_coor( IEN(ee, 1:n_en) );
    ux_ele = disp( IEN(ee, :),1 );
    uy_ele = disp( IEN(ee, :),2 );

    strain_zz=0;

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

    if question_def == 2 %plane stress
        stress=D0*[ strain_xx;strain_yy;2*strain_xy ];
        strain_zz=-v/E *( stress(1,1)+stress(2,1) );
        strain_zz_plot(ee)=strain_zz;
    end
    x_plot(ee)=x;y_plot(ee)=y;
    strain_xx_plot(ee)=strain_xx;strain_yy_plot(ee)=strain_yy;strain_xy_plot(ee)=strain_xy;
end
%云图代码
%xx
figure;
% 生成规则网格
[Xq, Yq] = meshgrid(linspace(min(x_plot), max(x_plot), n_el), linspace(min(x_plot), max(y_plot), n_el));

% 插值z到规则网格上
Zq = griddata(x_plot, y_plot,strain_xx_plot, Xq, Yq, 'v4');%v4虽然慢，但比cubic插值得到的结果更好

% 绘制云图
surf(Xq,Yq,Zq)
shading interp
colorbar;
%view([90, 90]);     % 调整试图位置
hold on;
title('Contour Plot of Function strain for xxh');
xlabel('x');
ylabel('y');

%yy
figure;
% 生成规则网格
[Xq, Yq] = meshgrid(linspace(min(x_plot), max(x_plot), n_el), linspace(min(x_plot), max(y_plot), n_el));

% 插值z到规则网格上
Zq = griddata(x_plot, y_plot,strain_yy_plot, Xq, Yq, 'v4');

% 绘制云图
surf(Xq,Yq,Zq)
shading interp
colorbar;
%view([90, 90]);     % 调整试图位置
hold on;
title('Contour Plot of Function strain for yyh');
xlabel('x');
ylabel('y');

%xy
figure;
% 生成规则网格
[Xq, Yq] = meshgrid(linspace(min(x_plot), max(x_plot), n_el), linspace(min(x_plot), max(y_plot), n_el));

% 插值z到规则网格上
Zq = griddata(x_plot, y_plot,strain_xy_plot, Xq, Yq, 'v4');

% 绘制云图
surf(Xq,Yq,Zq)
shading interp
colorbar;
%view([90, 90]);     % 调整试图位置
hold on;
title('Contour Plot of Function strain for xyh');
xlabel('x');
ylabel('y');

if question_def== 2 %zz
    figure;
    % 生成规则网格
    [Xq, Yq] = meshgrid(linspace(min(x_plot), max(x_plot), n_el), linspace(min(x_plot), max(y_plot), n_el));

    % 插值z到规则网格上
    Zq = griddata(x_plot, y_plot,strain_zz_plot, Xq, Yq, 'v4');

    % 绘制云图
    surf(Xq,Yq,Zq)
    shading interp
    colorbar;
    %view([90, 90]);     % 调整试图位置
    hold on;
    title('Contour Plot of Function strain for zzh');
    xlabel('x');
    ylabel('y');
end

%exact solution
%xx
figure;
x = linspace(min(x_plot), max(x_plot), 50);
y = linspace(min(x_plot), max(x_plot), 50);
[xq,yq]=meshgrid(x,y);
strain_xx=zeros(50,50); %和x,y所分份数有关
strain_xx(:,:)=f_strain_xx(xq,yq);

zq = griddata(x, y, strain_xx, xq, yq, 'cubic');
surf(xq,yq,zq)
shading interp
colorbar;
%view([90, 90]);     % 调整试图位置
hold on;
title('Contour Plot of Function strain for xx');
xlabel('x');
ylabel('y');

%yy
figure;
x = linspace(min(x_plot), max(x_plot), 50);
y = linspace(min(x_plot), max(x_plot), 50);
[xq,yq]=meshgrid(x,y);
strain_yy=zeros(50,50); %和x,y所分份数有关
strain_yy(:,:)=f_strain_yy(xq,yq);

zq = griddata(x, y, strain_yy, xq, yq, 'cubic');
surf(xq,yq,zq)
shading interp
colorbar;
%view([90, 90]);     % 调整试图位置
hold on;
title('Contour Plot of Function strain for yy');
xlabel('x');
ylabel('y');

%xy
figure;
x = linspace(min(x_plot), max(x_plot), 50);
y = linspace(min(x_plot), max(x_plot), 50);
[xq,yq]=meshgrid(x,y);
strain_xy=zeros(50,50); %和x,y所分份数有关
strain_xy(:,:)=f_strain_xy(xq,yq);

zq = griddata(x, y, strain_xy, xq, yq, 'cubic');
surf(xq,yq,zq)
shading interp
colorbar;
%view([90, 90]);     % 调整试图位置
hold on;
title('Contour Plot of Function strain for xy');
xlabel('x');
ylabel('y');

if question_def ==2 %%zz
    figure;
    x = linspace(min(x_plot), max(x_plot), 50);
    y = linspace(min(x_plot), max(x_plot), 50);
    [xq,yq]=meshgrid(x,y);
    strain_zz=zeros(50,50); %和x,y所分份数有关
    strain_zz(:,:)=f_strain_zz(xq,yq);

    zq = griddata(x, y, strain_zz, xq, yq, 'cubic');
    surf(xq,yq,zq)
    shading interp
    colorbar;
    %view([90, 90]);     % 调整试图位置
    hold on;
    title('Contour Plot of Function strain for zz');
    xlabel('x');
    ylabel('y');
end

