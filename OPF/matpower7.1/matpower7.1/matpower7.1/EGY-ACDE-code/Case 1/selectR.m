function [ Xr1 ,Xr2 ] = selectR( P,i )
%SELECTR 此处显示有关此函数的摘要
%   此处显示详细说明
    [m,~] = size(P);
    % PUA
    r1 = randi(m);
    while i==r1
        r1 = randi(m);
    end
    
    r2 = randi(m);
    while i==r2 || r1==r2
        r2 = randi(m);
    end
    
    Xr1 = P(r1,:);
    Xr2 = P(r2,:);   
end

