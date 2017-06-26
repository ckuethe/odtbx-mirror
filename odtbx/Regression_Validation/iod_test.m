function fail = iod_test
% Adaptor to run IODAppication_Tests in ODTBX regressionTesting harness
testCase = IODApplication_Tests;
Results  = run(testCase);
fail = any([Results.Failed]);