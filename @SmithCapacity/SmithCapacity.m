classdef SmithCapacity
    % SmithCapacity
    % A class of methods to compute the capacity according the Smith1971
    % and Smit1969.
    
    properties(Constant)
        % D is the noise entropy, which is 1/2 * log(2*pi*e)
        D = 0.5 * (log(2*pi)+1);
    end
    
    methods(Static)
        
        % Computes capacity according to the Smith1971 algorithm.
        [C, poi, voi] = computeCapacityOnlyAmplitudeCon(Alim, delta)
        [C, poi, voi] = computeCapacityOnlyAmplitudeConFast(Alim, nStart)
        
        [C, n] = computeTableOfCapacity(A)
        
        function r = H(poi, voi)
            % Computes the output entropy given F defined by the discrete
            % points of increase poi, with the corresponding values of
            % increase voi.  
            
            % Call xlogx instead of doing it in the anonymous function so
            % we don't need to compute p_Y twice for each value of y.  
            r = -integral( ...
                @(y) SmithCapacity.xlogx( ...
                SmithCapacity.p_Y(y, poi, voi)), -Inf, Inf, ...
                'AbsTol', 1e-12, ... % default AbsTol is 1e-10
                'RelTol', 1e-9); % default RelTol is 1e-6 (dominant)
        end % function H
        
        function r = I(poi, voi)
            % I(F) = H(F) - D
            r = SmithCapacity.H(poi, voi) - SmithCapacity.D;
        end % function I
        
        function r = i(x, poi, voi)
            % marginal information density
            f = @(x) integral( ...
                @(y) SmithCapacity.aloga_over_b( ...
                SmithCapacity.p_N(y-x), ...
                SmithCapacity.p_Y(y, poi, voi)), ...
                -Inf, Inf, 'AbsTol', 1e-12, 'RelTol', 1e-9);
            r = arrayfun(f, x);
        end % function i
        
        function r = I_Z(Z)
            % The information map I as a function of vector Z as defined in
            % Smith1971 p.212.
            %
            % Z is a vector of length 2*n.
            
            n = length(Z)/2;
            
            % break Z into voi and poi
            voi = Z(1:n);
            poi = Z(n+1:end);
            
            r = SmithCapacity.I(poi, voi);
        end % function I_Z
        
        optimal = checkCorollary1(A, poi, voi, I_Fo)
        
        function pdf = p_N(n)
            pdf = normpdf(n, 0, 1);
        end % function p_N
        
        function pdf = p_Y(y, poi, voi)
            % The pdf of Y = X + N, where X has points of increase at poi,
            % each with a voi increase in CDF.  
            % Parameter y may be a vector (either row or col).  
            % Output pdf is the same shape as y.
            
            % Store the original shape of y and unroll it into a column
            % vector.  
            orig_shape_y = size(y);
            y = y(:);
            
            % make poi a row vector
            if(iscolumn(poi))
                poi = poi';
            end
            % and voi a column vector
            if(isrow(voi))
                voi = voi';
            end
            
            M = bsxfun(@(y,p) normpdf(y, p, 1), y, poi);
            pdf = M * voi;
            
            % Restore the result to the shape of y
            pdf = reshape(pdf, orig_shape_y);
        end % function p_Y
        
    end % Methods(Static)
    
    % Internal methods (but publicly accessible to facilitate testing).
    methods(Static)
        function a = xlogx(x)
            a = x .* log(x);
            a(x==0) = 0;
        end % function xlogx

        function c = aloga_over_b(a, b)
            % computes a.*log(a./b)
            c = a.*log(a./b);
                        
            % Assume that if both a and b are zero, this is due to a
            % precision error (and that a and b are actually greater than
            % 0).  In this case, rather than return NaN, return the limit
            % of a*log(a/b), where
            % a = exp(-x^2)
            % b = exp(-k*x^2)
            % as x approaches 0.  According to Wolfram Alpha, this is 0.
            % Combined with c(a==0 & b>0) = 0,
            c(a==0) = 0;
            
            % Note that if(a>0 & b==0), then the result is NaN, and may
            % cause problems.  
            % This method is used in SmithCapacity.i to calculate i(x;F).
            % Examining the formula for i(x;F) in Smith1971 p.206, b should
            % only approach zero for extreme values of y, and in this case,
            % the limit as both of them approach zero (as explained above)
            % should be zero.  (Doing this to work around problems that
            % arise from precision errors.)
            c(b==0) = 0;
        end % function alog_a_over_b
    end % methods(Access = protected, Static)
    
end

