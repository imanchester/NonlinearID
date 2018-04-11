function dx = msd_nonlinearSpring_gen(x,u,t,param)
% Continuous time model of masses, connected by springs/dampers with
% params.

% New, more general version: handles different nonlinearities and model
% orders.

dx = zeros(size(x,1),1);

idt = param.idt;

% Control input:
uc = u(:,floor(t/idt)+1);

% To do: implement general nth order dynamics

% Hard-coded 4th order
% -------------------------------------------------------------------------

% Parameters:
m1 = param.m(1);
m2 = param.m(2);

c1 = param.c(1);
c2 = param.c(2);

k1 = param.k(1);
k2 = param.k(2);

nl1 = param.nl(1);
nl2 = param.nl(2);

slim1 = param.slim(1);
slim2 = param.slim(2);

slim_ = 1.25; % Recently changed from 10
smax = 1e5;

% Update velocity:
dx(1) = x(3);
dx(2) = x(4);

% Update acceleration:

% First compute the spring forces,
sf1 = nlspring(x(1),nl1,slim1);
sf2 = nlspring(x(1)-x(2),nl2,slim2); % Note: odd function

dx(3) = (-c1*x(3) - k1*sf1 - c2*(x(3)-x(4)) - k2*sf2 + uc(1))/(m1);
dx(4) = (k2*sf2 - c2*(x(4)-x(3)))/(m2);


%% Visualize the nonlinearity:
% figure
% r = -1.1:0.1:1.1;
% % r = -slim:0.1:slim;
% % plot(r,0.1*invsgmd(r))
% plot(r,nlspring(r))
% check = 1;

% Test a highly nonlinear spring defined by a bezier function
% figure
% r = 0:0.05:1.0;
% br = bezier(r,[0 1 0.5 0.0 1]);   
% plot(r,br)

% Test the final saturation behaviour
% figure
% r = 0:0.05:1.0;
% br = bezier(r,[1 2 5 1e3]);   
% plot(r,br)

% Now we need to combine these into a complete spring response:
% figure
% % r = -0.99*slim:0.05:0.99*slim;
% % r = 0:0.05:0.99*slim;
% r = 0:0.05:1.1;
% % r = 1.13:0.01:1.24;
% % r = -1.1:0.05:0;
% % r = -1.24:0.01:-1.13;
% sf = nlspring(r);   
% plot(r,sf)
% 
% dum = 1;


%% Nonlinear spring function

% Now accept different nonlinearity type argument and maximum displacement.

    function f = nlspring(a,nltype,slim)
       
        if a >= slim            
            f = smax;            
        elseif a <= -slim            
            f = -smax;            
        else      
  
% % Tan: better than inverse sigmoid
% % --------------------------------------------------------------------------
% % Bezier curve first:
%             s = abs(a)/slim;
%             b = sign(a).*(2*s.*(1-s)*1 + s.^2*(pi/2));            
% % Then pass through tan
%             f = tan(b);

% Tan: better than inverse sigmoid - used this in first promising results
% --------------------------------------------------------------------------
        if nltype == 1

% Bezier curve first:
            s = abs(a)/slim;
            b = sign(a).*(2*s.*(1-s)*1.5 + s.^2*(pi/2));            
% Then pass through tan
            f = tan(b);
            
% Completely custom Bezier: could be interesting!
% --------------------------------------------------------------------------                        
        elseif nltype == 2;
            
        str = 0.9;
        
        bc1 = [0 0 0.01 0.01 1 0.5 0.2 0.5 1]; bc2 = [1 2 5 1e3]; % Decent        
%         bc1 = [0 0.5 0.5 0.2 1]; bc2 = [1 2 5 1e3]; % Decent
%         bc1 = [0 1 0.5 0.0 1]; bc2 = [1 2 5 1e3]; % Not so much
%         bc1 = [0 2.5 0.5 -2.5 5]; bc2 = [5 10 20 1e3]; % Too extreme
        if abs(a) <= str*slim
% Apply initial force
            s = abs(a)/(str*slim); % Scale displacement
            f = sign(a).*bezier(s,bc1);
        else
% Apply saturation force
            s = (abs(a)-str*slim)/((1-str)*slim); % Scale displacement
            f = sign(a).*bezier(s,bc2);            
        end

        end % ...of nonlinearity types.
        
        end        
    end



    function s = invsgmd(a)
        
        if a >= slim            
            s = smax;            
        elseif a <= -slim            
            s = -smax;            
        else            
            s = -log((slim - a)./(a + slim));            
        end
        
    end


    function s = tannl(a)
        
        if a >= slim            
            s = smax;            
        elseif a <= -slim            
            s = -smax;            
        else            
            s = tan(a*pi/(2*slim));            
        end
    end

    function b = bezier(s,c)
       
% Inputs: 
%   - s: evaluation point in [0,1]
%   - c: vector of coefficients
    n = length(c) - 1;
    b = 0;
    for i = 0:n
        b = b + nchoosek(n,i)*((1-s).^(n-i)).*(s.^i)*c(i+1);
    end
        
    end

end

