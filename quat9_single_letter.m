% returns chaincodes, each row vector is a chaincode
close all;
clear all;
filename = 'test_single_letter_4';
output = quat9(filename);

function chaincodes = quat9(file)
    close all;
    data_sheet_path = 'C:\Users\nico\pen_proj\data\';
    folder = '';

    path = [data_sheet_path folder file '.csv'];
    chaincodes = process_file(path);
end

function chaincodes_ret = process_file(path)
    q_in = [];
    chaincodes = [];
    init_lens = [];
    new_lens = [];
    g_input = readtable(path);

    M = table2array(g_input);
    idx = all(isnan(M),2);
    idr = diff(find([1;diff(idx);1]));
    D = mat2cell(M,idr(:),size(M,2));
    D{1:2:end};

    j = 0;

    for m = 1:size(D)
        input_arr = D(m);
        in = cell2mat(input_arr);

        if size(in, 1) < 20 % TODO threshold, was 10
            continue;
        end
        j = j + 1;
%         if (j == 1) % skip first letter bc not part of recording
%             continue;
%         end
        
        quat_in = in(1:end,1:end) / 1073741824;
        q1 = quat_in(1:end,1);
        q2 = quat_in(1:end,2);
        q3 = quat_in(1:end,3);
        q0 = sqrt(abs(1.0 - ((q1.^2) + (q2.^2) + (q3.^2))));
        quat_in = [q0 q1 q2 q3];
        
%         print = ['num = ' num2str(j) '; len = ' num2str(size(q1,1))] % letter number and length
        [chaincode, init_len, new_len] = process_data(quat_in, j);
        init_lens = [init_lens; init_len];
        new_lens = [new_lens; new_len];
        
        valid_chaincode = init_len > 15; % at least 15 samples long
        if (valid_chaincode)
            chaincodes = [chaincodes; chaincode.'];
        else 
            j = j - 1;
        end
        
        q_in = [q_in; q0 q1 q2 q3];
%         break; % one letter only
    end
%     figure; hold on; grid on;
%     bar([1:1:size(init_lens,1)], init_lens);
%     hold off;
    
    chaincodes_ret = chaincodes;
    chaincodes_len = size(chaincodes, 1)
    lengths = [[1:1:size(init_lens,1)].', init_lens, new_lens]
    process_data(q_in, -1); % whole word
end

function [chaincode,init_len,new_len] = process_data(q, num)
    srt_i = 1;
    end_i = size(q(:,1),1);

    qw = q(:,1);
    qx = q(:,2);
    qy = q(:,3);
    qz = q(:,4);

    angle = 2 * acos(qw);
    x = qx ./ sqrt(1-qw.^2);
    y = qy ./ sqrt(1-qw.^2);
    z = qz ./ sqrt(1-qw.^2);
    
    x_i = x .* angle;
    y_i = y .* angle;
    z = z .* angle;

    y = x_i;
    x = y_i;
    
    % INITIAL SCALING
    a = 1;
    b = 5;
    [x, y] = scale_xy(x, y, a, b);
    z_scaled = (b - a) .* (z - min(z)) ./ (max(z) - min(z)) + a;
    z = z_scaled;
    
    % SMOOTHING
    degree = 5; % THRESHOLD moving avg
    x_fit = smooth(x, degree);
    y_fit = smooth(y, degree);
    z_fit = smooth(z, degree);
    
    % ROTATE Y-AXIS
    [x_fit,y_fit,z_fit] = rotate_y_axis(x_fit,y_fit,z_fit,-45); % angle in degrees

    x_fit = -x_fit;
    [x_fit, y_fit] = scale_xy(x_fit, y_fit, a, b);
    z_scaled = (b - a) .* (z_fit - min(z_fit)) ./ (max(z_fit) - min(z_fit)) + a;
    z_fit = z_scaled;
      
    init_len = size(x_fit,1);
    
%     figure; hold on; title('z over time'); grid on;
%     plot(z_fit, '-og');
%     hold off;
    
    figure; hold on; title(['letter: ' num2str(num) ', num points: ' num2str(init_len)]); grid on;
    subplot(2,2,1);
    plot(z_fit, '-og');
    subplot(2,2,2);
    plot(z_fit,x_fit, '-ob', 'LineWidth', 2);
    
    % truncate
    [z_fit,x_fit,zxgrad] = truncate_data(z_fit,x_fit,y_fit,num);
    new_len = size(z_fit,1);
    hold on;
    plot(z_fit,x_fit,'-og', 'LineWidth',2);
    
    P = polyfit(z_fit,x_fit,1);
    xlin = 1:1:5;
    ylin = P(1)*xlin+P(2);
    title(['fitted slope: ' num2str(P(1))])
    plot(xlin,ylin,'r-.');
    
    slope = P(1);
    if (slope < -0.5 && slope > -1)
        % rotate CW s.t. slope = vertical
        angle = atan(1/slope) * 180/pi / 2;
        print = ['letter: ' num2str(num) ', angle: ' num2str(angle)]
        [x_fit,y_fit,z_fit] = rotate_y_axis(x_fit,y_fit,z_fit,angle);
        [z_fit, x_fit] = scale_xy(z_fit, x_fit, a, b);
        plot(z_fit, x_fit, '-om', 'LineWidth', 2);
        P = polyfit(z_fit,x_fit,1);
        xlin = 1:1:5;
        ylin = P(1)*xlin+P(2);
        title(['fitted slope: ' num2str(P(1))])
        plot(xlin,ylin,'-r', 'LineWidth', 2);
    elseif (slope > 0.2)
        % rotate CCW
    end

%     hold off;
%     figure; title(['letter: ' num2str(num)]); hold on; grid on;
%     plot(zxgrad, '-ob');
%     hold off;
    
%     before_resampling = ['num = ' num2str(num), '; len = ' num2str(size(x_fit,1))]
    
    if (init_len < 15)
        chaincode = [];
        return;
    end
    
    [z_fit,x_fit] = normalize_char(z_fit,x_fit,num);
    
    zxgrad = gradient(x_fit) ./ gradient(z_fit);
%     figure; hold on; title(['letter: ' num2str(num)]); grid on;
    subplot(2,2,3);
    plot(zxgrad, '-ob');
%     hold off;
%     figure; hold on; title(['letter: ' num2str(num)]); grid on;
    subplot(2,2,4);
    h = histogram(zxgrad);
    hold off;
    
%     zx_dist = calc_dist(z_fit,x_fit)

    p1 = size(x_fit,1);
%     resample1 = ['num = ' num2str(num), '; len = ' num2str(p1)]
    % trajectory that will be used to derive chaincode
%     figure; hold on; grid on;
%     plot(z_fit,x_fit,'o-m', 'LineWidth', 2); 
%     hold off;
%     zx_dist = calc_dist(z_fit,x_fit) % distance isnt perfectly uniformly spaced but mostly uniform

    chaincode = chain_code(z_fit,x_fit);
    

end

function [ret_x,ret_y,xygrad] = truncate_data(x,y,z,num)
    start = x(1);
    if (start < 2.5) % new letter
        xgrad = gradient(x);
        cmp = xgrad < 0.02;
        temp = find(cmp);
        new_start = temp(1);
        x = x(new_start:end);
        y = y(new_start:end);
        [x, y] = scale_xy(x, y, 1, 5);
        
        start = x(1);
        if (start > x(end) && x(end) > 3)
            cmp = x < x(end);
            temp = find(cmp);
            if (size(temp,1) > 0)
                new_start = temp(1);
%                 print = num
                x = x(new_start:end);
                y = y(new_start:end);
                [x, y] = scale_xy(x, y, 1, 5);
            end
        end
        
        xgrad = gradient(x);
        ygrad = gradient(y);
        xcmp = xgrad < 0;
        ycmp = ygrad < 0;
        xtemp = find(~xcmp);
        ytemp = find(~ycmp);
        temp = setxor(xtemp,ytemp);
        if (size(temp,1) > 0)
            new_start = temp(1);
%             print2 = num
            x = x(new_start:end);
            y = y(new_start:end);
            [x, y] = scale_xy(x, y, 1, 5);
        end
        
        xygrad = ygrad ./ xgrad;
        
        ret_x = x;
        ret_y = y;
        return;
    end
    % new word   
    
    xlen = size(x,1);
    th = 60;
    if(xlen > th)
        x = x(xlen-th:end);
        y = y(xlen-th:end);
        [x, y] = scale_xy(x, y, 1, 5);
    end
    
    xgrad = gradient(x);
    cmp = xgrad > -0.02;
    temp = find(cmp);
    new_start = temp(1);
    x = x(new_start:end);
    y = y(new_start:end);
    [x, y] = scale_xy(x, y, 1, 5);
    
    print = ['letter: ' num2str(num) ', new length: ' num2str(size(x,1))];
    
    xgrad = gradient(x);
    ygrad = gradient(y);
    xygrad = ygrad ./ xgrad;

    ret_x = x;
    ret_y = y;
end

function chaincode = chain_code(x,y,num) % parameterizable by # of directions
    num_dir = 8; % number of directions, assumes divisible by 4 
    slice_angle = 2*pi / num_dir;

    x_dist = x(2:end) - x(1:end-1);
    y_dist = y(2:end) - y(1:end-1);
    angle = atan(x_dist./y_dist);
    chaincode = zeros(size(angle,1),1)-1;
    
    % upper half
    cmp = y_dist >= 0;
    temp = find(cmp);
    chaincode(temp) = mod(floor(angle(temp) / slice_angle),8);
    
    % lower half
    cmp = y_dist < 0;
    temp = find(cmp);
    chaincode(temp) = mod(4 + floor(angle(temp) / slice_angle),8);
    
%     size(chaincode)
end

function [ret_x,ret_y] = normalize_char(x,y,num)
    resample_len = 30; % THRESHOLD
    
    xy_dist = calc_dist(x,y);
    pt_dist = sum(xy_dist) / (resample_len);

    [z_norm,x_norm] = resample_char(x, y, pt_dist, resample_len);
%     plot(z_norm,x_norm, '-om', 'LineWidth', 1);
    
    z_temp = spline(1:1:size(z_norm,1), z_norm);
    x_temp = spline(1:1:size(x_norm,1), x_norm);
    xx = 0:size(z_norm,1)/(resample_len):size(z_norm,1);

    z_rs = ppval(z_temp,xx);
    z_rs = z_rs(2:end).';
    x_rs = ppval(x_temp,xx);
    x_rs = x_rs(2:end).';
    
    ret_x = z_rs;
    ret_y = x_rs;
    
end

function xy_dist = calc_dist(x,y)
    x_dist = abs(x(2:end) - x(1:end-1));
    y_dist = abs(y(2:end) - y(1:end-1));

    xy_dist = sqrt(x_dist.^2 + y_dist.^2);
end 

% some sort of smoothing, getting rid of extra points (assume less distance means less important) 
function [ret_x, ret_y] = resample_char(x,y,S,I)
    x = flip(x);
    y = flip(y);

    D = 0;
    p1 = [x(1),y(1)];
    ret_p = [p1];
    for i = 2:1:(size(x,1))
        p2 = [x(i),y(i)];
        dist = sqrt((p2(1) - p1(1))^2 + (p2(2) - p1(2))^2);
        num_seg = dist / S;
        if (num_seg < 1)
            continue;
        end
        if (num_seg == 1)
            ret_p = [ret_p; p2]; % append
            p1 = p2;
            continue;
        end
        num_seg = floor(num_seg);
        for j = 1:1:size(num_seg)
            dist = sqrt((p2(1) - p1(1))^2 + (p2(2) - p1(2))^2);
            % linear interpolate w p2
            p_eff(1) = p1(1) + S/dist * (p2(1) - p1(1));
            p_eff(2) = p1(2) + S/dist * (p2(2) - p1(2));
            ret_p = [ret_p; p_eff]; % append
            p1 = p_eff;
        end
    end % for each point in x
    ret_x = ret_p(:,1);
    ret_y = ret_p(:,2);
    
    ret_x = flip(ret_x);
    ret_y = flip(ret_y);
end

function [x_scaled,y_scaled] = scale_xy(x,y,a,b)
    x_scaled = (b - a) .* (x - min(x)) ./ (max(x) - min(x)) + a;
    y_scaled = (b - a) .* (y - min(y)) ./ (max(y) - min(y)) + a;
end

function [ret_x,ret_y,ret_z] = rotate_y_axis(x_fit,y_fit,z_fit,angle)
    % ROTATION y-AXIS
    R = roty(angle);
    Ar_mat = [];
    for i = 1:1:size(x_fit,1)
        Ab(1,1) = x_fit(i);
        Ab(2,1) = y_fit(i);
        Ab(3,1) = z_fit(i);
        Ar = (R*Ab).';
        Ar_mat = [Ar_mat; Ar];  
    end
    
    ret_x = Ar_mat(:,1);
    ret_y = Ar_mat(:,2);
    ret_z = Ar_mat(:,3);
end

function R = rotx(degrees)
    a = degrees * pi/180;
    R = [1,0,0 ; 0,cos(a),-sin(a) ; 0,sin(a),cos(a)];
end

function R = rotz(degrees)
    a = degrees * pi/180;
    R = [cos(a),-sin(a),0 ; sin(a),cos(a),0 ; 0,0,1];
end

function R = roty(degrees)
    a = degrees * pi/180;
    R = [cos(a),0,sin(a) ; 0,1,0 ; -sin(a),0,cos(a)];
end
