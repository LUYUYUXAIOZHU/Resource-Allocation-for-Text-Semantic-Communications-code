function [Matching,Cost] = Hungarian(Perf)
% 
% [MATCHING,COST] = Hungarian_New(WEIGHTS)
%
% 使用匈牙利算法在给定的MxN边权重矩阵WEIGHTS中寻找最小边权重匹配的函数
%
% 边权重为Inf表示该位置对应的顶点对没有相邻边
%
% MATCHING返回一个MxN矩阵，匹配位置为1，其他位置为0
% 
% COST返回最小匹配的成本

% 作者：Alex Melin，2006年6月30日


 % 初始化变量
 Matching = zeros(size(Perf));

% 通过移除任何未连接的顶点来压缩性能矩阵，以提高算法速度

  % 找出每列中连接的数量
    num_y = sum(~isinf(Perf),1);
  % 找出每行中连接的数量
    num_x = sum(~isinf(Perf),2);
    
  % 找出孤立的列（顶点）和行（顶点）
    x_con = find(num_x~=0);
    y_con = find(num_y~=0);
    
  % 组装压缩后的性能矩阵
    P_size = max(length(x_con),length(y_con));
    P_cond = zeros(P_size);
    P_cond(1:length(x_con),1:length(y_con)) = Perf(x_con,y_con);
    if isempty(P_cond)
      Cost = 0;
      return
    end

    % 确保存在完美匹配
      % 计算边矩阵的一种形式
      Edge = P_cond;
      Edge(P_cond~=Inf) = 0;
      % 找出边矩阵中的缺陷(CNUM)
      cnum = min_line_cover(Edge);
    
      % 投影额外的顶点和边，使得完美匹配存在
      Pmax = max(max(P_cond(P_cond~=Inf)));
      P_size = length(P_cond)+cnum;
      P_cond = ones(P_size)*Pmax;
      P_cond(1:length(x_con),1:length(y_con)) = Perf(x_con,y_con);
   
%*************************************************
% 主程序：控制执行哪个步骤
%*************************************************
  exit_flag = 1;
  stepnum = 1;
  while exit_flag
    switch stepnum
      case 1
        [P_cond,stepnum] = step1(P_cond);
      case 2
        [r_cov,c_cov,M,stepnum] = step2(P_cond);
      case 3
        [c_cov,stepnum] = step3(M,P_size);
      case 4
        [M,r_cov,c_cov,Z_r,Z_c,stepnum] = step4(P_cond,r_cov,c_cov,M);
      case 5
        [M,r_cov,c_cov,stepnum] = step5(M,Z_r,Z_c,r_cov,c_cov);
      case 6
        [P_cond,stepnum] = step6(P_cond,r_cov,c_cov);
      case 7
        exit_flag = 0;
    end
  end

% 移除所有虚拟的卫星和目标，并将匹配解压缩到原始性能矩阵的大小
Matching(x_con,y_con) = M(1:length(x_con),1:length(y_con));
Cost = sum(sum(Perf(Matching==1)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   步骤1：找出每行中最小的数字，并从该行中减去这个最小值
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [P_cond,stepnum] = step1(P_cond)

  P_size = length(P_cond);
  
  % 遍历每一行
  for ii = 1:P_size
    rmin = min(P_cond(ii,:));
    P_cond(ii,:) = P_cond(ii,:)-rmin;
  end

  stepnum = 2;
  
%**************************************************************************  
%   步骤2：在P_cond中找到一个零。如果其列或行中没有星标零，则标记该零。
%           对每个零重复此过程
%**************************************************************************

function [r_cov,c_cov,M,stepnum] = step2(P_cond)

% 定义变量
  P_size = length(P_cond);
  r_cov = zeros(P_size,1);  % 显示行是否被覆盖的向量
  c_cov = zeros(P_size,1);  % 显示列是否被覆盖的向量
  M = zeros(P_size);        % 显示位置是星标还是primed的掩码
  
  for ii = 1:P_size
    for jj = 1:P_size
      if P_cond(ii,jj) == 0 && r_cov(ii) == 0 && c_cov(jj) == 0
        M(ii,jj) = 1;
        r_cov(ii) = 1;
        c_cov(jj) = 1;
      end
    end
  end
  
% 重新初始化覆盖向量
  r_cov = zeros(P_size,1);  % 显示行是否被覆盖的向量
  c_cov = zeros(P_size,1);  % 显示列是否被覆盖的向量
  stepnum = 3;
  
%**************************************************************************
%   步骤3：用星标零覆盖每一列。如果所有列都被覆盖，则匹配是最大匹配
%**************************************************************************

function [c_cov,stepnum] = step3(M,P_size)

  c_cov = sum(M,1);
  if sum(c_cov) == P_size
    stepnum = 7;
  else
    stepnum = 4;
  end
  
%**************************************************************************
%   步骤4：找到一个未覆盖的零并prime它。如果包含这个primed零的行中没有星标零，
%           转到步骤5。否则，覆盖该行并取消包含星标零的列的覆盖。
%           继续这个过程直到没有未覆盖的零剩下。保存最小的未覆盖值并转到步骤6。
%**************************************************************************
function [M,r_cov,c_cov,Z_r,Z_c,stepnum] = step4(P_cond,r_cov,c_cov,M)

P_size = length(P_cond);

zflag = 1;
while zflag  
    % 找到第一个未覆盖的零
      row = 0; col = 0; exit_flag = 1;
      ii = 1; jj = 1;
      while exit_flag
          if P_cond(ii,jj) == 0 && r_cov(ii) == 0 && c_cov(jj) == 0
            row = ii;
            col = jj;
            exit_flag = 0;
          end      
          jj = jj + 1;      
          if jj > P_size; jj = 1; ii = ii+1; end      
          if ii > P_size; exit_flag = 0; end      
      end

    % 如果没有未覆盖的零，转到步骤6
      if row == 0
        stepnum = 6;
        zflag = 0;
        Z_r = 0;
        Z_c = 0;
      else
        % Prime未覆盖的零
        M(row,col) = 2;
        % 如果该行中有星标零
        % 覆盖该行并取消包含零的列的覆盖
          if sum(find(M(row,:)==1)) ~= 0
            r_cov(row) = 1;
            zcol = find(M(row,:)==1);
            c_cov(zcol) = 0;
          else
            stepnum = 5;
            zflag = 0;
            Z_r = row;
            Z_c = col;
          end            
      end
end
  
%**************************************************************************
% 步骤5：构建一系列交替的primed和星标零，如下所示。
%        让Z0表示在步骤4中找到的未覆盖primed零。
%        让Z1表示Z0列中的星标零（如果有）。
%        让Z2表示Z1行中的primed零（总是会有一个）。
%         继续直到系列在一个没有星标零在其列中的primed零处终止。
%         取消系列中每个星标零的星标，为系列中每个primed零加星标，
%         擦除所有primes并取消矩阵中每条线的覆盖。返回步骤3。
%**************************************************************************

function [M,r_cov,c_cov,stepnum] = step5(M,Z_r,Z_c,r_cov,c_cov)

  zflag = 1;
  ii = 1;
  while zflag 
    % 找到列中星标零的索引号
    rindex = find(M(:,Z_c(ii))==1);
    if rindex > 0
      % 保存星标零
      ii = ii+1;
      % 保存星标零的行
      Z_r(ii,1) = rindex;
      % 星标零的列与primed零的列相同
      Z_c(ii,1) = Z_c(ii-1);
    else
      zflag = 0;
    end
    
    % 如果primed零的列中有星标零，继续
    if zflag == 1;
      % 在最后一个星标零的行中找到primed零的列
      cindex = find(M(Z_r(ii),:)==2);
      ii = ii+1;
      Z_r(ii,1) = Z_r(ii-1);
      Z_c(ii,1) = cindex;    
    end    
  end
  
  % 取消路径中所有星标零的星标，并为所有primed零加星标
  for ii = 1:length(Z_r)
    if M(Z_r(ii),Z_c(ii)) == 1
      M(Z_r(ii),Z_c(ii)) = 0;
    else
      M(Z_r(ii),Z_c(ii)) = 1;
    end
  end
  
  % 清除覆盖
  r_cov = r_cov.*0;
  c_cov = c_cov.*0;
  
  % 移除所有primes
  M(M==2) = 0;

stepnum = 3;

% *************************************************************************
% 步骤6：将最小未覆盖值添加到每个覆盖行的每个元素，并从每个未覆盖列的每个元素中减去它。
%         返回步骤4，不改变任何星标、primes或覆盖线。
%**************************************************************************

function [P_cond,stepnum] = step6(P_cond,r_cov,c_cov)
a = find(r_cov == 0);
b = find(c_cov == 0);
minval = min(min(P_cond(a,b)));

P_cond(find(r_cov == 1),:) = P_cond(find(r_cov == 1),:) + minval;
P_cond(:,find(c_cov == 0)) = P_cond(:,find(c_cov == 0)) - minval;

stepnum = 4;

function cnum = min_line_cover(Edge)

  % 步骤2
    [r_cov,c_cov,M,stepnum] = step2(Edge);
  % 步骤3
    [c_cov,stepnum] = step3(M,length(Edge));
  % 步骤4
    [M,r_cov,c_cov,Z_r,Z_c,stepnum] = step4(Edge,r_cov,c_cov,M);
  % 计算缺陷
    cnum = length(Edge)-sum(r_cov)-sum(c_cov);