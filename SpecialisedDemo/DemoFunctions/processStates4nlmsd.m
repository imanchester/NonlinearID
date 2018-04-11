function [xd,yd,ud] = processStates4nlmsd(xsc,ysc,ud,dt)

% Apply a low pass filter

    % Filter frequency:
    fw = 0.1; %rad/s
    aw = 0.01; fw_rads = (1-aw)/(aw*dt); fw_hz = fw_rads/(2*pi);
    vel = dhpf(xsc(1:2,:),aw);
    vel = diag(max(abs(vel),[],2))\vel;

    xd = xsc;
    yd = ysc(:,1:size(xd,2));
    ud = ud(:,1:size(xd,2));

% Apply central differences 

    afb = 0.76;
    xf = dlpf(xsc(1:2,:),afb);
    xf = fliplr(dlpf(fliplr(xf),afb));

    vel = zeros(2,length(xd));

    % Central difference
    for t = 2:length(xf)-1
        vel(:,t) = xf(1:2,t+1) - xf(1:2,t-1);
    end
    vel(:,end) = xf(1:2,end) - xf(1:2,end-1);
    vel = diag(max(abs(vel),[],2))\vel;

    xd = [xsc(1:2,:);vel]; 


end

