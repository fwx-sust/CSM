function [n_np,n_el,x_coor,y_coor,IEN,top_pos,bottom_pos,left_pos,right_pos] = try_to_read_mesh(filename)
% 打开文件
run(filename);

fieldnames(msh);
POS=getfield(msh,'POS');
n_np=getfield(msh,'nbNod'); %节点个数
QUADS=getfield(msh,'QUADS');
n_el=size(QUADS,1); %number of element

%x_coor,y_coor
x_coor=POS(:,1,:);
y_coor=POS(:,2,:);

%IEN——QUADS
IEN = zeros(n_el, 4); %四边形element
for i=1:n_el
    IEN(i,1)=QUADS(i,1);
    IEN(i,2)=QUADS(i,2);
    IEN(i,3)=QUADS(i,3);
    IEN(i,4)=QUADS(i,4);
end

%boundary position
%5 for left, 6 for top, 7 for right, 8 for bottom
BC=getfield(msh,'LINES');
a=size(BC,1); %边界点个数
top_pos=[];
bottom_pos=[];
left_pos=[];
right_pos=[];

%right
flag=0;
for i=1:a
    if BC(i,3)==7
        flag=flag+1;
        right_pos(flag,1)=BC(i,1);
        right_pos(flag,2)=BC(i,2);
    end
end

%left
flag=0;
for i=1:a
    if BC(i,3)==5
        flag=flag+1;
        left_pos(flag,1)=BC(i,1);
        left_pos(flag,2)=BC(i,2);
    end
end

%bottom
flag=0;
for i=1:a
    if BC(i,3)==8
        flag=flag+1;
        bottom_pos(flag,1)=BC(i,1);
        bottom_pos(flag,2)=BC(i,2);
    end
end

%top
flag=0;
for i=1:a
    if BC(i,3)==6
        flag=flag+1;
        top_pos(flag,1)=BC(i,1);
        top_pos(flag,2)=BC(i,2);
    end
end