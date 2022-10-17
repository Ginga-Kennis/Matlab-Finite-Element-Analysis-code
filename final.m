%時間計測スタート
start_time = tic;
%vtkファイルを少し変えて読み込みやすいテキストファイルにし、読み込み
M = importdata('C:\Users\81801\Desktop\計算機3\Programing final\readgnesh3.txt');
num_nodes= M(1,1)             %接点数
num_elements = M(1,2)         %要素数
nodes = zeros(num_nodes,3);   %接点と要素行列の初期化
elements = zeros(num_elements,3);
for i = 2:num_nodes+1         %接点行列作成
    nodes(i-1,1) = M(i,1);
    nodes(i-1,2) = M(i,2);
    nodes(i-1,3) = M(i,3);
end
for i = 1:num_elements        %要素行列作成
    elements(i,1)= M(num_nodes+1+i,1);
    elements(i,2)= M(num_nodes+1+i,2);
    elements(i,3)= M(num_nodes+1+i,3);
end
bc_nodes = [1,2,3,4];       %境界条件の設定
bc_values = [0,100,0,100];
nodes
elements

%Ke,K,T,右辺ベクトルの初期化
ke = zeros(3,3);                
k = zeros(num_nodes,num_nodes);
t_array = zeros(num_nodes,1);
output_array = zeros(num_nodes,1);


%各三角形の係数行列の計算
for p = 1:num_elements
        v1 = elements(p,1)+1;
        v2 = elements(p,2)+1;
        v3 = elements(p,3)+1;

        x1 =nodes(v1,1);
        x2 =nodes(v2,1);
        x3 =nodes(v3,1);
        y1 =nodes(v1,2);
        y2 =nodes(v2,2);
        y3 =nodes(v3,2);

        a1=y2-y3;
        a2=y3-y1;
        a3=y1-y2;
        b1=x3-x2;
        b2=x1-x3;
        b3=x2-x1;
        s=(x1*y2+x2*y3+x3*y1-x1*y3-x2*y1-x3*y2)/2;
        %Keベクトルの作成
        ke = [a1*a1 + b1*b1 a1*a2+b1*b2 a1*a3+b1*b3;a2*a1+b2*b1 a2*a2+b2*b2 a2*a3+b2*b3;a3*a1+b3*b1 a3*a2+b3*b2 a3*a3+b3*b3]/(4*s);
        %Keから各値をとりKに代入
        for i = 1:3;
            for j = 1:3;
                k(elements(p,i)+1,elements(p,j)+1) = k(elements(p,i)+1,elements(p,j)+1) + ke(i,j);
            end
        end
end
%境界条件に合わせるためにKベクトルに0,1を追加
for i= 1:num_bcs;
    output_array(bc_nodes(i)+1) = bc_values(i);
    for p = 1:num_nodes;
        if p==bc_nodes(i)+1
            k(p,p) = 1;
        else 
            k(bc_nodes(i)+1,p) = 0;
        end
    end
end
%T（温度）の計算
t_array = inv(k)*output_array
%時間計測終了
end_time = toc(start_time)
disp(end_time)

