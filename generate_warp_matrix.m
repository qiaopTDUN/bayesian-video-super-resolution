function [Fw] = generate_warp_matrix(vx, vy, hsz)
    vx = reshape(vx,hsz);
    vy = reshape(vy,hsz);
    [gx, gy] = meshgrid(1:hsz(2), 1:hsz(1));
    gx = gx + vx;
    gy = gy + vy;
    
    gx = round(gx);
    gy = round(gy);
    gx(gx<1) = 1;
    gx(gx>hsz(2)) = hsz(2);
    gy(gy<1) = 1;
    gy(gy>hsz(1)) = hsz(1);
    
%     mx = gx<1 || gx>hsz(2);
%     my = gy<1 || gy>hsz(1);
    id = sub2ind(hsz,gy,gx);
    np = prod(hsz);
    Fw = sparse(1:np, id, ones(1,np), np, np); 
end