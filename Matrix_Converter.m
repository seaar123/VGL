function [Wvm,DimV]=Matrix_Converter(Wmv,W,Typ)
% Typ=1, W from M to V else, W from V to M
if Typ == 1
    [row,col,Wvm] = find(Wmv);
    [DimV xx]=size(Wvm);
else
    Wvm=W;
    indx_W=find(Wvm);
    [DimV xx]=size(indx_W);
    Wvm(indx_W)=Wmv;
end


