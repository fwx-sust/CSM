%位移可视化
x_plot=zeros(n_np,1);y_plot=zeros(n_np,1);x_disp_plot=zeros(n_np,1);y_disp_plot=zeros(n_np,1);
for i=1:n_np
    x_plot(i)=x_coor(i);
    y_plot(i)=y_coor(i);
    x_disp_plot(i)=disp(i,1);
    y_disp_plot(i)=disp(i,2);
end
% 生成规则网格
[Xq, Yq] = meshgrid(linspace(min(x_plot), max(x_plot), n_np), linspace(min(x_plot), max(y_plot), n_np));

% 插值z到规则网格上
Zq = griddata(x_plot, y_plot, x_disp_plot, Xq, Yq, 'cubic');%三次插值

% 绘制云图——uhx
figure;
contourf(Xq, Yq, Zq, 20, 'LineColor', 'none');
colorbar; % 显示颜色条
hold on;
title('Contour Plot of Function uhx');
xlabel('x');
ylabel('y');

%uhy
figure;
Zq= griddata(x_plot, y_plot, y_disp_plot, Xq, Yq, 'cubic');
contourf(Xq, Yq, Zq, 20, 'LineColor', 'none');
colorbar; % 显示颜色条
hold on;
title('Contour Plot of Function uhy');
xlabel('x');
ylabel('y');

%exact solution
%uy
figure;
x = linspace(min(x_plot), max(x_plot), 50);
y = linspace(min(x_plot), max(x_plot), 50);
[xq,yq]=meshgrid(x,y);
uy=zeros(50,50); %和x,y所分份数有关
uy(:,:)=fuy(xq,yq);

zq = griddata(x, y, uy, xq, yq, 'cubic');
contourf(xq, yq, zq, 20, 'LineColor', 'none');
colorbar; % 显示颜色条
hold on;
title('Contour Plot of Function uy');
xlabel('x');
ylabel('y');

%ux
figure;
ux=zeros(50,50);
ux(:,:)=fux(xq,yq);
zq = griddata(x, y, ux, xq, yq, 'cubic');
contourf(xq, yq, zq, 20, 'LineColor', 'none');
colorbar; % 显示颜色条
hold on;
title('Contour Plot of Function ux');
xlabel('x');
ylabel('y');