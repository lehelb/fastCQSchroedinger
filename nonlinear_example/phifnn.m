function [out,outx]=phifnn(a,lamf,x)
    out=greenfun(lamf,x+a)+greenfun(lamf,x-a);
    outx = -sqrt(lamf)*out;
end
