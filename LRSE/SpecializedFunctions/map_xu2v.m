function m = map_xu2v(x,u)

    nx = size(x,1);
    nu = size(u,1);
    
    m = zeros(nx+nu);
    
%   Install x elements 
    xind = 2:2:nx+nu+2;
    uind = 1:2:nx+nu+2;

    for i = 1:nx
        
        if i > nu
            m(nu+i,i) = 1;
        else
            m(xind(i),i) = 1;
        end
    end
    
%     for i = 1:nu
%         m(uind(i),nx+i) = 1;
%     end
    for i = 1:nu
        
        if i > nx
           m(nx+i,nx+i) = 1; 
        else
            m(uind(i),nx+i) = 1;
        end
    end

end