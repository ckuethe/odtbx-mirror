classdef rrdotang_tests < matlab.unittest.TestCase
    %Regression tests for angle-rate capability in rrdotang, an ODTBX
    %function
    
    
    properties
    end
    
    methods (Test)
        function testAzimuthRate(testCase)
            x1 = [0;0;0;0;0;0];
            x2 = [1;0;1;0;1;1];
            Y = x1 - x2;
            t = [0];
            
            % REQUIRED NUMJAC CODE
%             Fty = AzimuthDot(t,Y);
%             F = @AzimuthDot;
%             expSolution = numjac(F,t,Y,Fty,x1);
%             expSolution = expSolution(1,:);
%             A = isnan(expSolution);
%             for i = 1:length(A)
%                 if A(i)
%                     expSolution(i) = 0;
%                 end
%             end
            
            % Solution verified with numjac
            expy = 1;
            expH = [1 0 0 0 -1 0];
            options = odtbxOptions('measurement');
            options = setOdtbxOptions(options,'useRange',false);
            options = setOdtbxOptions(options,'useRangeRate',false);
            options = setOdtbxOptions(options,'useAngleRates',true);
            [acty, actH] = rrdotang(t,x1,x2,options);
            acty = acty(1);
            actH = actH(1,:);
            testCase.verifyEqual(acty, expy, 'AbsTol', 1e-20);
            testCase.verifyEqual(actH, expH, 'AbsTol', 1e-20);
        end
        function testElevationRate(testCase)
            x1 = [0;0;0;0;0;0];
            x2 = [1;1;0;0;0;1];
            Y = x1 - x2;
            t = [0];
            
            % REQUIRED NUMJAC CODE
%             Fty = ElevationDot(t,Y);
%             F = @ElevationDot;
%             expSolution = numjac(F,t,Y,Fty,x1);
%             expSolution = expSolution(1,:);
%             A = isnan(expSolution);
%             for i = 1:length(A)
%                 if A(i)
%                     expSolution(i) = 0;
%                 end
%             end
            
            % Solution verified with numjac
            expy = -0.707106781186547;
            expH = [-0.353553390593274 -0.353553390593274  0 0 0 0.707106781186547];
            options = odtbxOptions('measurement');
            options = setOdtbxOptions(options,'useRange',false);
            options = setOdtbxOptions(options,'useRangeRate',false);
            options = setOdtbxOptions(options,'useAngleRates',true);
            [acty, actH] = rrdotang(t,x1,x2,options);
            acty = acty(2);
            actH = actH(2,:);
            testCase.verifyEqual(acty, expy, 'AbsTol', 1e-14);
            testCase.verifyEqual(actH, expH, 'AbsTol', 1e-14);
        end
    end
    
end

