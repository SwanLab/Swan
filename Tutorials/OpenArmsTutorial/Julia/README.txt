
The main code is OpenArmsTutorial.jl. Some packages (around 10) need to be installed in Julia to be able to run the code. Don't worry, the REPL will tell you which packages to install, one by one. 

--------------------------------

PlotterNN: uses getOutput() multiple times from costFunction but it's not defined. 
	Option 1) Define a getOutput() function in 		  costNN
	Option 2) See if it's just the name that is  		  outdated and then substitude 		  		  it by the name of the correct		  		  function.