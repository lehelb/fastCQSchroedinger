function [out,outx]=phienn(a,lame,x)
    out=greenfun(lame,x+a)-greenfun(lame,x-a);
    outx = -sqrt(lame)*out;
end
    