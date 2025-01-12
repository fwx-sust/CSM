## 计算固体力学project

### gmsh工作（两个geo文件，四个.m文件）
  1. 通过gmsh生成了一个1*1矩形(可以通过进入.geo文件改变点的坐标来改变形状)和一个切除了1/4圆的矩形（可以通过进入.geo文件改变边长和半径来改变形状），并且划分了四边形网格（网格可以改变边界上的节点数），生成.geo文件，然后转换成为.m文件，再编写方法来读取这个.m文件。
  2. 对于完整矩形的读取，使用try_to_read_mesh.m的代码。对于切除了1/4圆的矩形，使用read_mesh_with_hole.m的代码。其中边界位置的读取依靠判断msh.LINES中第三列的值，所以如果Physical Curve的编号发生变化，需要在这两个代码中进行改变。nodal coordinates和IED的读取对于任何形状都适用。

### 代码部分
  1. drive.m: 对于1*1矩形（因为没设置circle边界）。可以在代码通过改变question_def来决定是平面应变问题还是平面应力问题；改变boundary_conditions矩阵来决定四条边的边界条件（Dirichlet or Neumann）；改变disp_top_f等公式改变Dirichlet条件，改变stress_top_f等公式改变Neumann条件。 ————其中test.m是manufactured solution，条件是plane strain问题，四条边x,y方向位移为0，指定位移理论解。
  2. vis_disp_for_square.m:位移可视化，所有图形均可，只要有d和理论解解析式（需要drive）
  3. vis_strain_for_square.m：应变可视化，（直角坐标系）所有图形均可，只要有d和理论解解析式（需要drive）
  4. vis_stress_for_square.m：应力可视化，（直角坐标系）所有图形均可，只要有d和理论解解析式（需要drive）
  5. error_for_square.m:计算error和他的convergence rates（无需drive，网格大小为代码定义）
  6. error_calculator.m：计算error，（直角坐标系）所有图形均可。（需要drive）
  7. stiffnessD.m：提供D矩阵
  9. drive_for_hole_Fig1: 带圆孔的Fig1 的主程序，但没写完，或许把body force和边界stress的函数补全了就对了呢QAQ
  10. drive_for_hole_Fig2：带圆孔的Fig2的主程序，同上QAQ，或许把body force补全了就对了呢QAQ
  11. vis_stress_for_hole：应力可视化，同square的不同是因为有坐标变换，所以理论解画图不太一样，但有限元解代码是一样的。受主程序的影响。
  12. 代码运行：运行drive，再运行其他想看的程序，但是drive里的限制和test中一样，只是可以自定义条件。对于manufactured solution，运行test.m再运行其他程序。
  
