function [ Xr1 ,Xr2 ] = selectR( P,i )
%SELECTR �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
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

