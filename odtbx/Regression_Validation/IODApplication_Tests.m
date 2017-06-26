classdef IODApplication_Tests < matlab.unittest.TestCase
    % IODAPPLICATION_TESTS Regression tests for IOD application
    %
    %   Contained within this test file are tests for each of the major
    %   functions included in the IOD application. They can be run
    %   individually or all at the same time. These tests should be run
    %   after any major changes are made to the source code for the IOD
    %   application. 
    %                To run all tests: 
    %                testCase = IODApplication_Tests
    %                Results  = run(testCase)
    %
    %                To run only one test:
    %                testCase = IODApplication_Tests
    %                Results  = run(testCase, 'testGibbs')
    %
    % REVISION HISTORY:
    %    Author          Date (MM/DD/YYYY)         Comment
    %    Ryan Willmot    08/04/2015                Original
    
    properties
        OriginalPath
    end
    
%     methods (TestMethodSetup)
%         function addSolverToPath(testCase)
%             testCase.OriginalPath = path; % assume running rom Regression_Validation dir
%             addpath(fullfile(fileparts(pwd),'..',filesep,'ODTBX_Source',filesep,'IOD'))
%         end
%         
%     end
%     
%     methods (TestMethodTeardown)
%         function restorePath(testCase)
%             path(testCase.OriginalPath)
%         end
%     end
    
    methods (Test)
        function testDouble_r(testCase)
            % Test case from Der (GEO) 
            % True position and velocity values
            % r2 = [40325.23 12355.73 -153.49]  (km)
            % v2 = [-0.8999 2.9391 0.0036]  (km/s)
            
            % NOTE: All Double-r advanced settings left at default values for this test:
            % Initial range estimate = 6378.137
            % Range step-size = 6378.137
            % Qmin = 150
            % Iterations = 500
            % Reduction factor = 2
            % Search Levels = 20
            
            global Type
            global radius
            global muglobal
            muglobal = 398600.4418;
            radius = 6378.137;
            Type = 'Earth';
            t1 = [2006, 09, 11, 04, 45, 44.073];
            t2 = [2006, 09, 11, 04, 50, 44.073];
            t3 = [2006, 09, 11, 04, 55, 44.073];
            time = [t1;t2;t3];
            r1 =[4993.01295987661; -2297.52999241947; 3225.16771430566];
            r2 = [5042.07579589666; -2187.76007642425; 3225.16771430566];
            r3 = [5088.72572520592; -2076.94318102115; 3225.16771430566];
            rsite = [r1 r2 r3];
            meas = [17.473272 -5.246774; 18.725895 -5.245155; 19.978504 -5.243423];
            fprintf(1,'\n');
            [actr2, actv2] = Double_r(meas, time, rsite);
            expr2 = [40704.883405069; 9901.40160953186; -231.732898612314];
            expv2 = [-0.718225598133399; 2.96765960752106; 0.00338902557395866];
            testCase.verifyEqual(actr2, expr2, 'AbsTol', 1e-9);
            testCase.verifyEqual(actv2, expv2, 'AbsTol', 1e-13);
            clearvars -global;    
        end
        function testDouble_r_Accuracy(testCase)
            M32 = 0.174469646069451;
            M21 = 0.404537300198097;
            n = 0.000443237629661882;
            t1 = [2012,08,20,11,40,28];
            t2 = [2012,08,20,11,48,28];
            t3 = [2012,08,20,11,52,28];
            time = [t1;t2;t3];
            actSolution = Double_r_Accuracy(M32,M21,n,time);
            expSolution = 459.150296262386;
            testCase.verifyEqual(actSolution, expSolution, 'AbsTol', 1e-11);
        end
        function testDouble_r_OrbitElements(testCase)
            global Type
            global radius
            global muglobal
            muglobal = 398600.4418;
            radius = 6378.137;
            Type = 'Earth';
            t1 = [2006, 09, 11, 04, 45, 44.073];
            t2 = [2006, 09, 11, 04, 50, 44.073];
            t3 = [2006, 09, 11, 04, 55, 44.073];
            time = [t1;t2;t3];
            r1 =[4993.01295987661; -2297.52999241947; 3225.16771430566];
            r2 = [5042.07579589666; -2187.76007642425; 3225.16771430566];
            r3 = [5088.72572520592; -2076.94318102115; 3225.16771430566];
            rsite = [r1 r2 r3];
            R1mag = 37895.1030351563;
            R2mag = 37895.1030351563;
            LOS1 = [.949860536690582; 0.299002796225035; -0.0914455503974318];
            LOS2 = [0.943099596101956; 0.319696744804357; -0.0914174118744659];
            LOS3 = [0.935888121004651; 0.340237835432138; -0.0913873093095426];
            LOS = [LOS1 LOS2 LOS3];
            [~, M32, M21, n, R2, R3, E32, v2] = Double_r_OrbitElements(R1mag,R2mag,rsite,LOS,time);
            expM32 = 0.0184750807683453;
            expM21 = 0.0184715857341072;
            expn = 7.52264941296388e-05;
            expR2 = [36903.1443761328; 8612.66858793995; 136.780823675319];
            expR3 = [36707.3343312763; 9417.85700936602; 137.683751459251];
            expE32 = 0.020133937226519;
            expv2 = 0;
            testCase.verifyEqual(M32,expM32, 'AbsTol', 1e-15);
            testCase.verifyEqual(M21,expM21, 'AbsTol', 1e-15);
            testCase.verifyEqual(n,expn, 'AbsTol', 1e-18);
            testCase.verifyEqual(R2,expR2, 'AbsTol', 1e-9);
            testCase.verifyEqual(R3,expR3, 'AbsTol', 1e-10);
            testCase.verifyEqual(E32,expE32, 'AbsTol', 1e-15);
            testCase.verifyEqual(v2,expv2, 'AbsTol', eps);
            clearvars -global;
        end
        function testDouble_r_RangeEstimates(testCase)
            global Type
            global radius
            global muglobal
            muglobal = 398600.4418;
            radius = 6378.137;
            Type = 'Earth';
            t1 = [2006, 09, 11, 04, 45, 44.073];
            t2 = [2006, 09, 11, 04, 50, 44.073];
            t3 = [2006, 09, 11, 04, 55, 44.073];
            time = [t1;t2;t3];
            r1 =[4993.01295987661; -2297.52999241947; 3225.16771430566];
            r2 = [5042.07579589666; -2187.76007642425; 3225.16771430566];
            r3 = [5088.72572520592; -2076.94318102115; 3225.16771430566];
            rsite = [r1 r2 r3];
            meas = [17.473272 -5.246774; 18.725895 -5.245155; 19.978504 -5.243423];
            [actr1o, actr2o] = Double_r_RangeEstimates(meas, time, rsite);
            expr1o = 37895.1030351563;
            expr2o = 37895.1030351563;
            testCase.verifyEqual(actr1o,expr1o,'AbsTol',1e-9);
            testCase.verifyEqual(actr2o,expr2o,'AbsTol',1e-9);
            clearvars -global
        end
        function testGauss(testCase)
            % Test case from Vallado, 4th Edition page 447
            % True position and velocity values
            % r2 = [6356.486 5290.532 6511.396]  (km)
            % v2 = [-4.1729 4.7766 1.7203]  (km/s)
            
            global Type
            global radius
            global muglobal
            muglobal = 398600.4418;
            radius = 6378.137;
            Type = 'Earth';
            t1 = [2012,08,20,11,40,28];
            t2 = [2012,08,20,11,48,28];
            t3 = [2012,08,20,11,52,28];
            time = [t1;t2;t3];
            r1 = [4054.881; 2748.195; 4074.237];
            r2 = [3956.224; 2888.232; 4074.364];
            r3 = [3905.073; 2956.935; 4074.430];
            rsite = [r1 r2 r3];
            meas = [0.939913 18.667717; 45.025748 35.664741; 67.886655 36.996583];
            fprintf(1,'\n');
            [actr2, actv2] = Gauss(meas, time, rsite);
            expr2 = [6313.37813100503; 5247.50563424429; 6467.70716523848];
            expv2 = [-4.18548828180492; 4.78849291837557; 1.72171466021237];
            testCase.verifyEqual(actr2, expr2, 'AbsTol', 1e-10);
            testCase.verifyEqual(actv2, expv2, 'AbsTol', 1e-13);
            clearvars -global
        end
        function testGeneralLagrange(testCase)
            t1 = [2012,08,20,11,40,28];
            t2 = [2012,08,20,11,48,28];
            t3 = [2012,08,20,11,52,28];
            time = [t1;t2;t3];
            r1 = [4054.881; 2748.195; 4074.237];
            r2 = [3956.224; 2888.232; 4074.364];
            r3 = [3905.073; 2956.935; 4074.430];
            r = [r1 r2 r3];
            actSolution = GeneralLagrange(time, r);
            rexp = [3956.224; 2888.232; 4074.364];
            rdotexp = [-0.210597916666667;0.288089583333333;0.000271527777777199];
            r2dotexp = [-2.10937500000052e-05;-1.52256944444436e-05;2.89351851831851e-08];
            expSolution = [rexp rdotexp r2dotexp];
            testCase.verifyEqual(actSolution, expSolution, 'AbsTol', 1e-14);
        end
        function testGibbs(testCase)
            global muglobal
            muglobal = 398600.4418;
            r1 =[8004.7213; 2812.996; 5408.8835 ];
            r2 = [6313.3958; 5247.5237; 6467.7250];
            r3 = [5272.0417; 6321.1254; 6810.4754];
            r = [r1 r2 r3];
            [actSolutionv2, actSolutione] = Gibbs(r);
            expSolutionv2 = [-4.18551770539569;4.78853260973447;1.72173226422994];
            expSolutione = 0.197847212874192;
            testCase.verifyEqual(actSolutionv2, expSolutionv2, 'AbsTol', 1e-13);
            testCase.verifyEqual(actSolutione, expSolutione, 'AbsTol', 1e-14);
            clearvars -global
        end
        function testHerrick_Gibbs(testCase)
            global muglobal
            muglobal = 398600.4418;
            t1 = [2012,08,20,11,40,28];
            t2 = [2012,08,20,11,48,28];
            t3 = [2012,08,20,11,52,28];
            time = [t1;t2;t3];
            r1 =[8004.7213; 2812.996; 5408.8835 ];
            r2 = [6313.3958; 5247.5237; 6467.7250];
            r3 = [5272.0417; 6321.1254; 6810.4754];
            r = [r1 r2 r3];
            actSolution = Herrick_Gibbs(r, time);
            expSolution = [-4.10497457412767; 4.69686758382983; 1.68895090670568];
            testCase.verifyEqual(actSolution, expSolution, 'AbsTol', 1e-13);
            clearvars -global
        end
        function testJDate(testCase)
            t1 = [2012,08,20,11,40,28];
            t2 = [2012,08,20,11,48,28];
            t3 = [2012,08,20,11,52,28];
            time = [t1;t2;t3];
            actSolution = JDate(time);
            expSolution = [2456159.98643519 2456159.99199074 2456159.99476852];
            testCase.verifyEqual(actSolution, expSolution, 'RelTol', 1e-13);
        end
        function testLaplace(testCase)
            % Test case from Vallado, 4th Edition page 447
            % True position and velocity values
            % r2 = [6356.486 5290.532 6511.396]  (km)
            % v2 = [-4.1729 4.7766 1.7203]  (km/s)
            
            global Type
            global radius
            global muglobal
            muglobal = 398600.4418;
            radius = 6378.137;
            Type = 'Earth';
            t1 = [2012,08,20,11,40,28];
            t2 = [2012,08,20,11,48,28];
            t3 = [2012,08,20,11,52,28];
            time = [t1;t2;t3];
            r1 = [4054.881; 2748.195; 4074.237];
            r2 = [3956.224; 2888.232; 4074.364];
            r3 = [3905.073; 2956.935; 4074.430];
            rsite = [r1 r2 r3];
            meas = [0.939913 18.667717; 45.025748 35.664741; 67.886655 36.996583];
            fprintf(1,'\n');
            [actr2, actv2] = Laplace(meas, time, rsite);
            expr2 = [6405.71357461414; 5339.92410393526; 6561.46022107998];
            expv2 = [-4.2944952378819; 4.16491065422315; 1.26796787461234];       
            testCase.verifyEqual(actr2, expr2, 'AbsTol', 1e-10);
            testCase.verifyEqual(actv2, expv2, 'AbsTol', 1e-13);
            clearvars -global
        end
        function testLOS_Vectors(testCase)
            meas = [0.939913 18.667717; 45.025748 35.664741; 67.886655 36.996583];
            actSolution = LOS_Vectors(meas);
            LOS1 = [0.9472633016838; 0.0155408474202159; 0.320079239165193];
            LOS2 = [0.57422536042227; 0.574741691745486; 0.583041356352574];
            LOS3 = [0.300651905274714; 0.739921912374617; 0.601767393136729];
            expSolution = [LOS1 LOS2 LOS3];
            testCase.verifyEqual(actSolution, expSolution, 'AbsTol', 1e-14);
        end
        function testLST(testCase)
            date = [1996 08 20 8 30 0];
            long = -110;
            actSolution = LST(date, long);
            expSolution = 346.456265351425;
            testCase.verifyEqual(actSolution, expSolution, 'AbsTol', 1e-11);
        end
        function testsigpt(testCase)
            global Type
            global muglobal
            global radius
            Type = 'Sun';
            muglobal = 132712400180.6;
            radius = 695800;
            t1 = [2014, 05, 24, 0, 0, 0];
            t2 = [2014, 05, 25, 0, 0, 0];
            t3 = [2014, 05, 26, 0, 0, 0];
            time = [t1; t2; t3];
            x = [247.373905; 2.834073; 247.485666; 2.613774; 247.596334; 2.389592];
            f = @Gauss;
            r1 = [-0.4647054599604308;-0.9005772369452507;-8.277565505977721e-05];
            r2 = [-0.44964265600355746;-0.9084425422680832;-8.206114809721906e-05];
            r3 = [-0.43444965806586244;-0.9160490874486091;-8.132840124312375e-05];
            rsite = [r1 r2 r3] * 149597871;
            Px1 = [3.2828064058084e-09;0;0;0;0;0];
            Px2 = [0;3.2828064058084e-09;0;0;0;0];
            Px3 = [0;0;3.2828064058084e-09;0;0;0];
            Px4 = [0;0;0;3.2828064058084e-09;0;0];
            Px5 = [0;0;0;0;3.2828064058084e-09;0];
            Px6 = [0;0;0;0;0;3.2828064058084e-09];
            Px = [Px1 Px2 Px3 Px4 Px5 Px6];
            [acty, actPy] = sigpt(x,Px,f,time,rsite);
            expy = [-85167555.2869454;-179089679.084418;2121974.14549082;28.3069677563881;-11.0404517687539;-2.23857661982802];
            Py1 = [3523960781073.69 8501408518761.59 -420083451030.772 -414828.4073396 -460986.829918809 440743.081093577];
            Py2 = [8501408518761.59 20509293751504.1 -1013433816860.16 -1000753.60444168 -1112105.14678391 1063274.00055244];
            Py3 = [-420083451030.772 -1013433816860.16 50077206603.9823 49452.6492264268 54957.755721683 -52540.2176419765];
            Py4 = [-414828.4073396 -1000753.60444168 49452.6492264268 0.052248639130557 0.0625078650775927 -0.0522895058233086];
            Py5 = [-460986.829918809 -1112105.14678391 54957.755721683 0.0625078650775927 0.0801886426529095 -0.0586371322445382];
            Py6 = [440743.081093577 1063274.00055244 -52540.2176419765 -0.0522895058233086 -0.0586371322445382 0.0551724768007207];
            expPy = [Py1;Py2;Py3;Py4;Py5;Py6];
            % Tolerance values are larger for this test case because output
            % values are very large
            testCase.verifyEqual(acty, expy, 'RelTol', 1e-13);
            testCase.verifyEqual(actPy, expPy, 'RelTol', 1e-13);
            clearvars -global
        end
        function testSite_Position(testCase)
           global radius
           global Type
           radius = 6378.137;
           Type = 'Earth';
           LST =  346.456265351425;
           h_ELLP = 2;
           lat = 40;
           actSolution = Site_Position(lat, LST, h_ELLP);
           expSolution = [4758.13796088895;-1146.16985121307;4079.27113108534];
           testCase.verifyEqual(actSolution, expSolution, 'AbsTol', 1e-10);
           clearvars -global
        end
    end
    
end

