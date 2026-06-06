classdef B_SphHarmonic < handle
    % BphHarmonic - B 向量球面谐波，对应于wiki百科中的第2个函数,
    % 对应于kargl JASA 1993 附录中的B函数,
    % https://en.wikipedia.org/wiki/Vector_spherical_harmonics#cite_note-1
    % 对应于Morse & Feshbach《Methods of theoretical physics Part.II》p.1865中的B函数。
    properties (SetAccess='private')
        parent; % Parent object
        % Bx,By,Bz在一个立方体上的B向量
        Bx;
        By;
        Bz;
        % 在球面上
        Bx_sph;
        By_sph;
        Bz_sph;
        Bln_sph;
        % Bav:由curl()产生的B的角速度
        Bav;
        Bangx;
        Bangy;
        Bangz;
        Bang;   
    end
    
    methods
        %% 1.设置函数
        function hObj = B_SphHarmonic(parent)
            hObj.parent = parent;
            % 计算矢量球谐函数 B = curl(C)/k
            % [hObj.Bx, hObj.By, hObj.Bz, hObj.Bav] = ...
            %     curl(hObj.parent.x, hObj.parent.y, hObj.parent.z,...
            %     hObj.parent.C.Cx, hObj.parent.C.Cy, hObj.parent.C.Cz);
            % % 这里应该除以 k，但由于我们使用的是球谐函数而不是适当的函数，所以我们只需将波长设置为 1。
            % hObj.Bx = hObj.Bx/(2*pi);
            % hObj.By = hObj.By/(2*pi);
            % hObj.Bz = hObj.Bz/(2*pi);

            % 计算矢量球谐函数 B = r×C;
            hObj.Bx = hObj.parent.y.*hObj.parent.C.Cz   - hObj.parent.z.*hObj.parent.C.Cy;      %（r(vector)×▽Yln）在x方向分量
            hObj.By = -(hObj.parent.x.*hObj.parent.C.Cz - hObj.parent.z.*hObj.parent.C.Cx);     % ...         y..
            hObj.Bz = hObj.parent.x.*hObj.parent.C.Cy   - hObj.parent.y.*hObj.parent.C.Cx;

            % 获取单位球面上的矢量分量
            hObj.Bx_sph = interp3(hObj.parent.x, hObj.parent.y, hObj.parent.z, hObj.Bx,...
                                  hObj.parent.x_sph, hObj.parent.y_sph, hObj.parent.z_sph);
            hObj.By_sph = interp3(hObj.parent.x, hObj.parent.y, hObj.parent.z, hObj.By,...
                                  hObj.parent.x_sph, hObj.parent.y_sph, hObj.parent.z_sph);
            hObj.Bz_sph = interp3(hObj.parent.x, hObj.parent.y, hObj.parent.z, hObj.Bz,...
                                  hObj.parent.x_sph, hObj.parent.y_sph, hObj.parent.z_sph);
            hObj.Bln_sph = sqrt( hObj.Bx_sph.^2 + hObj.By_sph.^2 + hObj.Bz_sph.^2 );

            % 需要将B与单位径向向量取叉积
            hObj.Bangx = hObj.parent.y.*hObj.Bz   - hObj.parent.z.*hObj.By;
            hObj.Bangy = -(hObj.parent.x.*hObj.Bz - hObj.parent.z.*hObj.Bx);
            hObj.Bangz = hObj.parent.x.*hObj.By   - hObj.parent.y.*hObj.Bx;

            % 获取单位球面上的矢量分量
            hObj.Bangx = interp3(hObj.parent.x, hObj.parent.y, hObj.parent.z,...
                                 hObj.Bangx, hObj.parent.x_sph, hObj.parent.y_sph, hObj.parent.z_sph);
            hObj.Bangy = interp3(hObj.parent.x, hObj.parent.y, hObj.parent.z,...
                                 hObj.Bangy, hObj.parent.x_sph, hObj.parent.y_sph, hObj.parent.z_sph);
            hObj.Bangz = interp3(hObj.parent.x,hObj.parent.y, hObj.parent.z,...
                                 hObj.Bangz, hObj.parent.x_sph, hObj.parent.y_sph, hObj.parent.z_sph);
            hObj.Bang  = sqrt(hObj.Bangx.^2 + hObj.Bangy.^2 + hObj.Bangz.^2);
        
        end
        
        %% 2.绘图函数  
        function plot_ang(hObj)
            figure;
            surf(hObj.parent.x_sph.*real(hObj.Bang), hObj.parent.y_sph.*real(hObj.Bang), ...
                hObj.parent.z_sph.*real(hObj.Bang), real(hObj.Bang),...
                'EdgeColor', 'flat', 'FaceColor','interp');
            switch hObj.parent.parity 
                case Parity.Even
                    % title_str = sprintf('$$|\\mathrm{(r/|r|)\\times B^e_{%d,%d}}|$$', hObj.parent.n, hObj.parent.m);
                    title_str = sprintf('$$|\\mathrm{B_{%d,%d}^e}|$$', hObj.parent.n, hObj.parent.m);
                case Parity.Odd
                    % title_str = sprintf('$$|\\mathrm{(r/|r|)\\times B^o_{%d,%d}}|$$', hObj.parent.n, hObj.parent.m);
                    title_str = sprintf('$$|\\mathrm{B_{%d,%d}^o}|$$', hObj.parent.n, hObj.parent.m);
                otherwise
                    % title_str = sprintf('$$|\\mathrm{(r/|r|)\\times B_{%d,%d}}|$$'  , hObj.parent.n, hObj.parent.m);
                    title_str = sprintf('$$|\\mathrm{B_{%d,%d}}|$$'  , hObj.parent.n, hObj.parent.m);
            end
            title(title_str,'Interpreter','latex');
            xlabel('x');
            ylabel('y');
            zlabel('z');
            cb = colorbar;
            colormap jet
            ylabel(cb, title_str , 'interpreter' ,'latex');
        end
        
        function plot_ang_vec(hObj)
            figure;
            quiver3(hObj.parent.x_sph, hObj.parent.y_sph, hObj.parent.z_sph, ...
                hObj.Bangx, hObj.Bangy, hObj.Bangz, hObj.parent.arrow_scale,...
                'Color', 'black');
            switch hObj.parent.parity 
                case Parity.Even
                    title_str = sprintf('$$\\mathrm{(r/|r|)\\times B^e_{%d,%d}}$$', hObj.parent.n, hObj.parent.m);
                case Parity.Odd
                    title_str = sprintf('$$\\mathrm{(r/|r|)\\times B^o_{%d,%d}}$$', hObj.parent.n, hObj.parent.m);
                otherwise
                    title_str = sprintf('$$\\mathrm{(r/|r|)\\times B_{%d,%d}}$$', hObj.parent.n, hObj.parent.m);
            end
            title(title_str,'Interpreter','latex');
            xlabel('x');
            ylabel('y');
            zlabel('z');
            hold on;
            surf(hObj.parent.x_sph, hObj.parent.y_sph, hObj.parent.z_sph, real(hObj.parent.Yln_sph),...
                    'EdgeColor', 'interp', 'FaceColor','none');
            hold off;
            colormap jet;
            pbaspect([1 1 1]);
        end
        
        function plot_ang_abs_vec(hObj)
            figure;
            % 向量场是 B 的实部，因此我们想要向量实部的绝对值。
            abs_vec = sqrt(real(hObj.Bangx).^2 + real(hObj.Bangy).^2 + real(hObj.Bangz).^2);
            surf(hObj.parent.x_sph.*abs_vec, hObj.parent.y_sph.*abs_vec, hObj.parent.z_sph.*abs_vec,...
                abs_vec, 'EdgeColor', 'flat', 'FaceColor','interp');
            
            switch hObj.parent.parity 
                case Parity.Even
                    title_str = sprintf('$$\\mathrm{|(r/|r|)\\times B^e_{%d,%d}|}$$', hObj.parent.n, hObj.parent.m);
                case Parity.Odd
                    title_str = sprintf('$$\\mathrm{|(r/|r|)\\times B^o_{%d,%d}|}$$', hObj.parent.n, hObj.parent.m);
                otherwise
                    title_str = sprintf('$$\\mathrm{|(r/|r|)\\times B_{%d,%d}|}$$', hObj.parent.n, hObj.parent.m);
            end
            title(title_str,'Interpreter','latex');
            xlabel('x');
            ylabel('y');
            zlabel('z');
            hold on;
            surf(hObj.parent.x_sph, hObj.parent.y_sph, hObj.parent.z_sph, real(hObj.parent.Yln_sph),...
                    'EdgeColor', 'interp', 'FaceColor','none');
            hold off;
            colormap jet;
            pbaspect([1 1 1]);
        end

        %% 其它函数
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

