classdef C_SphHarmonic < handle
    % CphHarmonic - C 向量球面谐波，对应于wiki百科中的第3个矢量球谐函数
    % 对应于kargl JASA 1993 附录中的C函数,
    % https://en.wikipedia.org/wiki/Vector_spherical_harmonics#cite_note-1
    % 对应于Morse & Feshbach《Methods of theoretical physics Part.II》p.1865中的C函数
    properties (SetAccess='private')
        parent; 
        % Cx,Cy,Cz -> 空间立方体上的 C 矢量
        Cx;
        Cy;
        Cz;
        % 球体表面上的C矢量
        Cx_sph;
        Cy_sph;
        Cz_sph;
        Cln_sph;       
        Cav; % 由 curl() 产生的 C 的角速度
        Cangx;
        Cangy;
        Cangz; % 仅 C 的 phi 和 theta 分量。
        Cang   % 模值
    end
    
    methods
        %% 设置函数
        function hObj = C_SphHarmonic(parent)
            % Calculate the vector spherical harmonic C=curl(rPhi)
            % 计算矢量球谐函数 C = ▽×(rYln) = r×▽Yln
            hObj.parent = parent;
            % [hObj.Cx, hObj.Cy, hObj.Cz, hObj.Cav] = ...
            %     curl(hObj.parent.x, hObj.parent.y, hObj.parent.z, ...
            %     hObj.parent.x.*hObj.parent.Yln, hObj.parent.y.*hObj.parent.Yln, hObj.parent.z.*hObj.parent.Yln);
            [grad_x, grad_y, grad_z] = gradient(hObj.parent.Yln, 2/hObj.parent.num_steps); % ▽Yln
            hObj.Cx = hObj.parent.y.*grad_z   - hObj.parent.z.*grad_y;      %（r(vector)×▽Yln）在x方向分量
            hObj.Cy = -(hObj.parent.x.*grad_z - hObj.parent.z.*grad_x);     % ...         y..
            hObj.Cz = hObj.parent.x.*grad_y   - hObj.parent.y.*grad_x;      % ...         z..
            hObj.Cx = hObj.Cx / sqrt( hObj.parent.n * (hObj.parent.n + 1) );
            hObj.Cy = hObj.Cy / sqrt( hObj.parent.n * (hObj.parent.n + 1) );
            hObj.Cz = hObj.Cz / sqrt( hObj.parent.n * (hObj.parent.n + 1) );

            % 得到单位球面上的矢量分量
            hObj.Cx_sph = interp3(hObj.parent.x, hObj.parent.y, hObj.parent.z, ...
                                  hObj.Cx, hObj.parent.x_sph, hObj.parent.y_sph, hObj.parent.z_sph);
            hObj.Cy_sph = interp3(hObj.parent.x, hObj.parent.y, hObj.parent.z, ...
                                  hObj.Cy, hObj.parent.x_sph, hObj.parent.y_sph, hObj.parent.z_sph);
            hObj.Cz_sph = interp3(hObj.parent.x, hObj.parent.y, hObj.parent.z, ...
                                  hObj.Cz, hObj.parent.x_sph, hObj.parent.y_sph, hObj.parent.z_sph);
            hObj.Cln_sph = sqrt( hObj.Cx_sph.^2 + hObj.Cy_sph.^2 + hObj.Cz_sph.^2);          
            
            % 求C与单位径向矢量的叉乘
            hObj.Cangx = hObj.parent.y   .* hObj.Cz   - hObj.parent.z .* hObj.Cy;
            hObj.Cangy = -(hObj.parent.x .* hObj.Cz   - hObj.parent.z .* hObj.Cx);
            hObj.Cangz = hObj.parent.x   .* hObj.Cy   - hObj.parent.y .* hObj.Cx;
    
            % 得到单位球面上的矢量分量
            hObj.Cangx = interp3(hObj.parent.x, hObj.parent.y, hObj.parent.z,...
                                 hObj.Cangx, hObj.parent.x_sph, hObj.parent.y_sph, hObj.parent.z_sph,'spline');
            hObj.Cangy = interp3(hObj.parent.x, hObj.parent.y, hObj.parent.z,...
                                 hObj.Cangy, hObj.parent.x_sph, hObj.parent.y_sph, hObj.parent.z_sph,'spline');
            hObj.Cangz = interp3(hObj.parent. x,hObj.parent.y, hObj.parent.z,...
                                 hObj.Cangz, hObj.parent.x_sph, hObj.parent.y_sph, hObj.parent.z_sph,'spline');
            hObj.Cang = sqrt(hObj.Cangx.^2 + hObj.Cangy.^2 + hObj.Cangz.^2);
        end
        
        %% 画图函数
        function plot_ang(hObj)
            figure;
            surf(hObj.parent.x_sph.*real(hObj.Cang), hObj.parent.y_sph.*real(hObj.Cang), ...
                hObj.parent.z_sph.*real(hObj.Cang), real(hObj.Cang),...
                'EdgeColor', 'flat', 'FaceColor','interp');     
            % surf(hObj.parent.x_sph.*imag(hObj.Cang), hObj.parent.y_sph.*imag(hObj.Cang), ...
            %     hObj.parent.z_sph.*imag(hObj.Cang), real(hObj.Cang),...
            %     'EdgeColor', 'flat', 'FaceColor','interp');
            switch hObj.parent.parity 
                case Parity.Even
                    % title_str = sprintf('$$|\\mathrm{(r/|r|)\\times C^e_{%d,%d}}|$$', hObj.parent.n, hObj.parent.m);
                    title_str = sprintf('$$|\\mathrm{C_{%d,%d}^e}|$$', hObj.parent.n, hObj.parent.m);
                case Parity.Odd
                    % title_str = sprintf('$$|\\mathrm{(r/|r|)\\times C^o_{%d,%d}}|$$', hObj.parent.n, hObj.parent.m);
                    title_str = sprintf('$$|\\mathrm{C_{%d,%d}^o}|$$', hObj.parent.n, hObj.parent.m);
                otherwise
                    % title_str = sprintf('$$|\\mathrm{(r/|r|)\\times C_{%d,%d}}|$$', hObj.parent.n, hObj.parent.m);
                    title_str = sprintf('$$|\\mathrm{C_{%d,%d}}|$$', hObj.parent.n, hObj.parent.m);
            end
            title(title_str,'Interpreter','latex');
            xlabel('x');
            ylabel('y');
            zlabel('z');
            cb = colorbar;
            colormap jet;
            ylabel(cb, title_str, 'interpreter', 'latex');
        end
        
        function plot_ang_vec(hObj)
            figure
            quiver3(hObj.parent.x_sph, hObj.parent.y_sph, hObj.parent.z_sph, ...
                real(hObj.Cangx), real(hObj.Cangy), real(hObj.Cangz), hObj.parent.arrow_scale,...
                'Color', 'black', 'LineWidth', 1);
            pbaspect([1 1 1]);
            switch hObj.parent.parity 
                case Parity.Even
                    title_str = sprintf('$$\\mathrm{(r/|r|)\\times C^e_{%d,%d}}$$', hObj.parent.n, hObj.parent.m);
                case Parity.Odd
                    title_str = sprintf('$$\\mathrm{(r/|r|)\\times C^o_{%d,%d}}$$', hObj.parent.n, hObj.parent.m);
                otherwise
                    title_str = sprintf('$$\\mathrm{(r/|r|)\\times C_{%d,%d}}$$', hObj.parent.n, hObj.parent.m);
            end
            title(title_str,'Interpreter','latex');
            xlabel('x');
            ylabel('y');
            zlabel('z');
            hold on;
            %surf(hObj.parent.x_sph, hObj.parent.y_sph, hObj.parent.z_sph, real(hObj.parent.Yln_sph),...
            %        'EdgeColor', 'interp', 'FaceColor','none');
            hold off;
            colormap jet;
        end
        
        function plot_ang_abs_vec(hObj)
            figure;
            % 向量场是复数 C 的实部，因此我们想要向量实部的绝对值。
            abs_vec = sqrt(real(hObj.Cangx).^2 + real(hObj.Cangy).^2 + real(hObj.Cangz).^2);
            surf(hObj.parent.x_sph.*abs_vec, hObj.parent.y_sph.*abs_vec, hObj.parent.z_sph.*abs_vec,...
                abs_vec, 'EdgeColor', 'flat', 'FaceColor','interp');
            pbaspect([1 1 1]);
            %hObj.set_pbaspect2;
            switch hObj.parent.parity 
                case Parity.Even
                    title_str = sprintf('$$\\mathrm{|(r/|r|)\\times C^e_{%d,%d}|}$$', hObj.parent.n, hObj.parent.m);
                case Parity.Odd
                    title_str = sprintf('$$\\mathrm{|(r/|r|)\\times C^o_{%d,%d}|}$$', hObj.parent.n, hObj.parent.m);
                otherwise
                    title_str = sprintf('$$\\mathrm{|(r/|r|)\\times C_{%d,%d}|}$$', hObj.parent.n, hObj.parent.m);
            end
            title(title_str,'Interpreter','latex');
            xlabel('x');
            ylabel('y');
            zlabel('z');
            colormap jet;
        end
        
        % 其他函数
        function set_pbaspect(hObj)
            switch hObj.parent.n
                case 0
                    pba = [1 1 1];
                case 1
                    switch hObj.parent.m
                        case -1
                            pba = [2 1 1];
                        case 0
                            pba = [1 1 2];
                        case 1
                            pba = [2 1 1];
                    end
                case 2
                    switch hObj.parent.m
                        case -2
                            pba = [3 3 1];
                        case -1
                            pba = [3 1 3];
                        case 0
                            pba = [1 1 4];
                        case 1
                            pba = [3 1 3];
                        case 2
                            pba = [3 3 1];
                    end
                case 3
                    switch hObj.parent.m
                        case -3
                            pba = [4 4 1];
                        case -2
                            pba = [3 3 2];
                        case -1
                            pba = [3 1 4];
                        case 0
                            pba = [1 1 6];
                        case 1
                            pba = [3 1 4];
                        case 2
                            pba = [3 3 2];
                        case 3
                            pba = [4 4 1];
                    end
                otherwise
                    pba = [1 1 1];
            end
            if hObj.parent.n<=3
                pbaspect(pba);
            end
        end
    end
end

