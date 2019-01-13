function RunCoverageTest
fileName = '/home/alex/git-repos/FEM-MAT-OO/TestsInMatlabFormat.m';
folderCovPath = 'FEM/PostProcess/Printer/ResultsPrinter';
codeCov = matlab.unittest.plugins.CodeCoveragePlugin.forFolder(folderCovPath,'IncludingSubfolders',true);
testToRun   = matlab.unittest.TestSuite.fromFile(fileName);
runner  = matlab.unittest.TestRunner.withTextOutput;
runner.addPlugin(codeCov);
result = runner.run(testToRun)
end

