% returns chaincodes, each row vector is a chaincode
% close all;
% clear all;
% filename = 'ax_data3';
% output = quat9(filename);

% ccode = [0;1;2;3;4;5;6;0;7];
% ccode = [0;1;2;3;4;5;6;7;8;9;10;11;12;13;14;15];
% ccode = [12;13;14;15;0];
% ccode = [0;1;2;3;4;5;6;7;8;9;10;11;12;13;14;15;16;17;18;19];
% ccode = 0:1:23;
% chaincode2trajectory(ccode.',24)

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
    zero_cnts = [];
    g_input = readmatrix(path); %readtable(path);
    M = g_input;
    %M = table2array(g_input);
    idx = all(isnan(M),2);
    idr = diff(find([1;diff(idx);1]));
    D = mat2cell(M,idr(:),size(M,2));
    D{1:2:end};

    j = 0;

    for m = 1:size(D)
        input_arr = D(m);
        in = cell2mat(input_arr);

        if size(in, 1) < 20
            continue;
        end
        j = j + 1;
        
        quat_in = in(1:end,1:end) / 1073741824;
        q1 = quat_in(1:end,1);
        q2 = quat_in(1:end,2);
        q3 = quat_in(1:end,3);
        q0 = sqrt(abs(1.0 - ((q1.^2) + (q2.^2) + (q3.^2))));
        acc_x = in(1:end,4);
        quat_in = [q0 q1 q2 q3 acc_x];
        
%         print = ['num = ' num2str(j) '; len = ' num2str(size(q1,1))] % letter number and length
        [chaincode, init_len, zero_cnt] = process_data(quat_in, j);
        init_lens = [init_lens; init_len];
        zero_cnts = [zero_cnts; zero_cnt];
        
        valid_chaincode = init_len > 15; % at least 15 samples long
        if (valid_chaincode)
            chaincodes = [chaincodes; chaincode.'];
        else 
            j = j - 1;
        end
        
        q_in = [q_in; q0 q1 q2 q3 acc_x];
    end
    
    chaincodes_ret = chaincodes;
    chaincodes_len = size(chaincodes, 1)
    lengths = [[1:1:size(init_lens,1)].', init_lens, zero_cnts]
%     process_data(q_in, -1); % whole word
end

function [chaincode,init_len,new_len] = process_data(q, num)
    srt_i = 1;
    end_i = size(q(:,1),1);

    qw = q(:,1);
    qx = q(:,2);
    qy = q(:,3);
    qz = q(:,4);
    ax = q(:,5);

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
    ax = smooth(ax,degree);
    
    % ROTATE Y-AXIS
    [x_fit,y_fit,z_fit] = rotate_y_axis(x_fit,y_fit,z_fit,-45); % angle in degrees

    x_fit = -x_fit;
    [x_fit, y_fit] = scale_xy(x_fit, y_fit, a, b);
    z_scaled = (b - a) .* (z_fit - min(z_fit)) ./ (max(z_fit) - min(z_fit)) + a;
    z_fit = z_scaled;
      
    init_len = size(x_fit,1);
    
    len = init_len;
    ax_scaled = (b - a) .* (ax - min(ax)) ./ (max(ax) - min(ax)) + a;
    
    figure;
    subplot(1,2,1); hold on; grid on;
    title(['letter: ' num2str(num) ', num points: ' num2str(init_len)]);
    plot(z_fit, '-ob', 'LineWidth', 2);
    plot(x_fit, '-or', 'LineWidth', 2);
    plot(ax_scaled, '-og', 'LineWidth', 2);
%     hold off;
    subplot(1,2,2); hold on; grid on;
    plot(z_fit,x_fit, '-ob', 'LineWidth', 2);
    
    % truncate
    [z_fit,x_fit] = truncate_char(z_fit,x_fit,ax_scaled,num);
    new_len = size(z_fit,1);
    title(['letter: ' num2str(num) ', num points: ' num2str(init_len) ', new len: ' num2str(new_len)]);
    subplot(1,2,2); hold on; grid on;
    plot(z_fit,x_fit, '-og', 'LineWidth', 2);
    
    [z_fit,x_fit] = normalize_char(z_fit,x_fit);
    
    chaincode = chain_code(z_fit,x_fit);
    hold off;
end

function [ret_x,ret_y] = truncate_char(x,y,ax,num)
    [ax_max, ax_prom] = islocalmax(ax);
    cmp = ax_prom > 0;
    temp = find(cmp);
    subplot(1,2,1);
    yyaxis right;
    plot(temp, ax_prom(temp), '.m', 'MarkerSize', 20);
    
    cmp = ax_prom > 0.2; %.4
    temp = find(cmp);
    cmp2 = temp > 10;
    temp2 = find(cmp2);
    if (size(temp2,1) < 1)
       ret_x = x;
       ret_y = y;
       return;
    end
    start_sample = temp(temp2(1));
    % new word
    if(x(1) > 3.5 && size(temp2,1) > 1)
        start_sample = temp(temp2(2));
%         print1 = [num, size(x,1) - start_sample + 1]
        if ((size(x,1) - start_sample + 1) > 50 && size(temp2,1) > 2)
           start_sample = temp(temp2(3));
%            print2 = num
        end
    end
    
    yyaxis left;
    plot(start_sample, [x(start_sample) ax(start_sample)], 'xm', 'MarkerSize', 20, 'LineWidth', 2);
    
    % truncate and scale
    x_fit = x(start_sample:end); % truncate
    y_fit = y(start_sample:end); % truncate
    [ret_x, ret_y] = scale_xy(x_fit, y_fit, 1, 5);
end

function chaincode = chain_code(x,y) % parameterizable by # of directions
    num_dir = 32; % number of directions, assumes divisible by 4 
    slice_angle = 2*pi / num_dir;

    x_dist = x(2:end) - x(1:end-1);
    y_dist = y(2:end) - y(1:end-1);
    angle = atan(x_dist./y_dist);
    chaincode = zeros(size(angle,1),1)-1;
    
    % upper half
    cmp = y_dist >= 0;
    temp = find(cmp);
    chaincode(temp) = mod(floor(angle(temp) / slice_angle),num_dir);
    
    % lower half
    cmp = y_dist < 0;
    temp = find(cmp);
    chaincode(temp) = mod(num_dir/2 + floor(angle(temp) / slice_angle),num_dir);
    
    chaincode2trajectory(chaincode,num_dir);
    
%     size(chaincode)
end

function chaincode2trajectory(ccode,ndir)
    ndir_div4 = ndir/4;
    y1 = 1;
    x1 = 1;
    ploty = y1;
    plotx = x1;
    for num = ccode.'
        signy = 1;
        signx = 1;
        if (num >= ndir/2)
            signx = -1;
        end
        if ((num >= ndir_div4) && (num <= (3*ndir_div4)))
            signy = -1;
        end
        m = mod(num, ndir_div4);
        if (m == 0)
            if (num == 0 || num == 2*ndir_div4)
                da = 0;
                db = 1;
            elseif (num == ndir_div4 || num == 3*ndir_div4)
                da = 1;
                db = 0; 
            end
            y1 = y1 + db*signy;
            x1 = x1 + da*signx;
        else 
            tan_angle = 2*pi / ndir; % (pi/2) / (ndir/4)
            if (signy*signx < 0)
               m = ndir_div4 - m;
            end
            slope = tan(m*tan_angle);
        
            % solve for unsigned da and db
            da = sqrt(slope^2 / (1+slope^2));
            db = sqrt(1-da^2);
            if (signy*signx < 0)
                y1 = y1 + db*signy;
                x1 = x1 + da*signx;
                num;
            else
                y1 = y1 + db*signy;
                x1 = x1 + da*signx;
            end
        end
        ploty = [ploty; y1];
        plotx = [plotx; x1];
    end
    [plotx, ploty] = scale_xy(plotx, ploty, 1, 5);    
    plot(plotx, ploty, '-or', 'LineWidth', 2);
%     xy_dist = calc_dist(plotx,ploty)
end

function [ret_x,ret_y] = normalize_char(x,y,num)
    resample_len = 30; % THRESHOLD
    
    xy_dist = calc_dist(x,y);
    pt_dist = sum(xy_dist) / (resample_len);

    [z_norm,x_norm] = resample_char(x, y, pt_dist, resample_len);
    subplot(1,2,2);
    plot(z_norm,x_norm, '-om', 'LineWidth', 1);
    
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
