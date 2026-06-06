classdef A_SphHarmonic < handle
    % SphHarmonic - 矢量球面谐波函数,对应于wiki百科中的第一个函数,
    % https://en.wikipedia.org/wiki/Vector_spherical_harmonics#cite_note-1
    % 对应于kargl JASA 1993 附录中的A函数,
    % 对应于Morse & Feshbach《Methods of theoretical physics Part.II》p.1865中的A函数。
    properties (SetAccess = 'private')
        parity; % 谐波的奇偶性
        n;  % 球谐函数的次数
        m;  % 球谐函数的阶数
        Yln; % 在空间立方体区域上定义的标量球谐函数
        Yln_sph; % 在单位球面上定义的标量球谐函数
        C; % M 球谐函数
        B; % N 球谐函数
        % 空间立方体区域的 x,y,z 坐标。
        x;
        y;
        z;
        num_steps;
        % 三维空间中的Theta 和 Phi
        Theta;
        Phi;
        R;
        % 球体的笛卡尔坐标
        x_sph;
        y_sph;
        z_sph;
        % 球体上的 theta 和 phi 值
        theta_sph;
        phi_sph;
        % 绘图的矢量箭头长度
        arrow_scale = 3;
    end

    methods
        % 1.设置函数
        function hObj = A_SphHarmonic(radius, m, n, sigma, num_steps)
            % 输入参数：
            % radius: 障碍物半径
            % n：球谐函数的指标 n
            % m：球谐函数的指标 m (m=0,1,2,...n)
            % sigma: 奇偶指标
            % num_steps: 空间网格步长
            switch sigma
                case 'e'
                    parity = Parity.Even;
                case 'o'
                    parity = Parity.Odd;
                case 'c'
                    parity = Parity.Complex;
                otherwise
                    error('invalid parity: sigma must be one of even , odd or complex.');
            end

            if nargin == 4
                hObj.num_steps = 30;
            else
                hObj.num_steps = num_steps;
            end

            if n < 0 || m < 0 || m > n
                ME = MException('Illegal value for n or m');
                throw(ME)
            end

            % 赋值
            hObj.n = n;
            hObj.m = m;
            hObj.parity = parity;
            % 球体表面上的 theta 和 phi 值(r=1)
            hObj.theta_sph = 0 : pi/(2*hObj.num_steps) : pi;
            hObj.phi_sph   = 0 : pi/hObj.num_steps     : 2*pi;
            % 障碍物球体的球坐标参数方程
            hObj.x_sph = radius * sin(hObj.theta_sph) .* cos(hObj.phi_sph)';
            hObj.y_sph = radius * sin(hObj.theta_sph) .* sin(hObj.phi_sph)';
            hObj.z_sph = radius * cos(hObj.theta_sph) .* ones(size(hObj.phi_sph))';
            % 定义笛卡尔坐标系上从 -r 到 r 的方形空间区域,计算此区域的球谐函数值再插值到单位球面上
            [hObj.x, hObj.y, hObj.z] = meshgrid( ...
                -radius : 2/hObj.num_steps*radius : radius, ...
                -radius : 2/hObj.num_steps*radius : radius, ...
                -radius : 2/hObj.num_steps*radius : radius ...
                );
            % 三维空间中球坐标系下的 Theta, Phi, R
            [hObj.Phi, hObj.Theta, hObj.R] = cart2sph(hObj.x, hObj.y, hObj.z);

            hObj.createYlm();
            hObj.C = C_SphHarmonic(hObj); % 球谐函数C
            hObj.B = B_SphHarmonic(hObj); % 球谐函数B
        end

        function createYlm(hObj) % 创建标量球谐函数
            Lln = legendre( hObj.n, cos(hObj.Theta + pi/2) );
            % warning: cart2sph返回的Theta为(-pi/2,pi/2),球坐标系的俯仰角范围为(0,pi)。
            if hObj.n ~= 0
                Lln = squeeze( Lln(hObj.m+1,:,:,:) );
            end

            if hObj.m == 0 % εm
                epsilon = 1;
            else
                epsilon = 2;
            end

            a1 = ( epsilon * (2*hObj.n+1) / (4*pi) );
            a2 = factorial( hObj.n - hObj.m ) / factorial( hObj.n + hObj.m );
            Const  = sqrt(a1*a2);

            switch hObj.parity
                case Parity.Even
                    hObj.Yln = Const*Lln.*cos(hObj.m*hObj.Phi);
                case Parity.Odd
                    hObj.Yln = Const*Lln.*sin(hObj.m*hObj.Phi);
                otherwise
                    hObj.Yln = Const*Lln.*exp(-1i*hObj.m*hObj.Phi);
            end
            hObj.Yln_sph = interp3( ... % 标量球谐函数在单位球上的值(注意这里是3维插值2维)
                hObj.x, hObj.y, hObj.z, hObj.Yln, ...
                hObj.x_sph, hObj.y_sph, hObj.z_sph);
        end

        % 直接设置l, m, parity的函数
        function set_l_m_parity(hObj, n, m, parity)
            if n < 0 || m < 0 || m > n
                ME = MException('SphHarmonic: Illegal value for n or m');
                throw(ME)
            end
            hObj.n = n;
            hObj.m = m;
            hObj.parity = parity;
            hObj.createYlm();
            hObj.C = C_SphHarmonic(hObj);
            hObj.B = B_SphHarmonic(hObj);
        end

        % 2.绘图函数
        function plot_Ylm(hObj) % 绘制球谐函数的强度随空间坐标的变化图
            Ylns = interp3( hObj.x, hObj.y, hObj.z, hObj.Yln, hObj.x_sph, hObj.y_sph, hObj.z_sph );
            % [Xr,Yr,Zr] = sph2cart(Phi2,Theta2-pi/2,real(Ymns).^2);
            figure
            surf( hObj.x_sph.*real(Ylns).^2, hObj.y_sph.*real(Ylns).^2, hObj.z_sph.*real(Ylns).^2, real(Ylns),...
                'EdgeColor', 'flat', 'FaceColor','interp');
            % surf( hObj.x_sph, hObj.y_sph, hObj.z_sph, real(Ylns),...
            %         'EdgeColor', 'flat', 'FaceColor','interp');
            hObj.set_pbaspect();
            switch hObj.parity
                case Parity.Even
                    title_str = sprintf('$|Y_{%d%d}^e|$', hObj.n, hObj.m);
                case Parity.Odd
                    title_str = sprintf('$|Y_{%d%d}^o|$', hObj.n, hObj.m);
                otherwise
                    title_str = sprintf('Y_%d^%d', hObj.n, hObj.m);
            end
            cb = colorbar;
            ylabel(cb, title_str , 'interpreter' ,'latex');
            xlabel('x');
            ylabel('y');
            zlabel('z');
            title(title_str, 'interpreter' ,'latex');
            colormap jet;
            axis equal
        end

        function plot_Cang(hObj) % 绘制矢量球谐函数Cσmn的切向分量模值分布图
            hObj.C.plot_ang();
        end

        function plot_CangVec(hObj) % 绘制矢量球谐函数Cσmn的切向分量矢量图
            hObj.C.plot_ang_vec();
        end

        function plot_CRadiation(hObj) % C型矢量球谐函数的切向分量绝对值辐射模式图
            hObj.C.plot_ang_abs_vec();
        end

        function plot_Bang(hObj) % 绘制矢量球谐函数Bσmn的切向分量模值分布图
            hObj.B.plot_ang();
        end

        function plot_BangVec(hObj) % 绘制矢量球谐函数Bσmn的切向分量矢量图
            hObj.B.plot_ang_vec();
        end

        function plot_BRadiation(hObj) % B型矢量球谐函数的切向分量绝对值辐射模式图
            hObj.B.plot_ang_abs_vec();
        end

        % 其他函数
        function set_pbaspect(hObj)
            % 设置球谐函数图的轴，使得每个轴具有相同的比例。s.t.所有轴都不会拥挤。
            switch hObj.n
                case 0
                    pba = [1 1 1];
                case 1
                    switch hObj.m
                        case -1
                            pba = [2 1 1];
                        case 0
                            pba = [1 1 2];
                        case 1
                            pba = [2 1 1];
                    end
                case 2
                    switch hObj.m
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
                    switch hObj.m
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
            if hObj.n <= 3
                pbaspect(pba);
            end
        end
    end
end