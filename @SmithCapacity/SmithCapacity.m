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
        
        function r = H(poi, voi)
            % Computes the output entropy given F defined by the discrete
            % points of increase poi, with the corresponding values of
            % increase voi.  
            
            % Call xlogx instead of doing it in the anonymous function so
            % we don't need to compute p_Y twice for each value of y.  
            r = -integral( ...
                @(y) SmithCapacity.xlogx( ...
                SmithCapacity.p_Y(y, poi, voi)), -Inf, Inf);
        end % function H
        
        function r = I(poi, voi)
            % I(F) = H(F) - D
            r = SmithCapacity.H(poi, voi) - SmithCapacity.D;
        end % function I
        
        function r = i(x, poi, voi)
            r = zeros(size(x));
            for ii = 1:numel(x)
                % marginal information density
                r(ii) = integral( ...
                    @(y) SmithCapacity.aloga_over_b( ...
                    SmithCapacity.p_N(y-x(ii)), ...
                    SmithCapacity.p_Y(y, poi, voi)), ...
                    -Inf, Inf);
            end % for ii
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
            
            % Parameter y may be a vector (of unknown shape).  Rather than
            % try to neatly vectorize everything across each poi too, just
            % implement this in a for loop.
            
            pdf = zeros(size(y));
            for ii = 1:length(poi)
                pdf = pdf + voi(ii)*normpdf(y, poi(ii), 1);
            end % for ii
        end % function p_Y
        
    end % Methods(Static)
    
    methods(Access = protected, Static)
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
