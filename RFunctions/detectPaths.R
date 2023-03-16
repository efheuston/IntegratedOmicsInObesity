# Setting up a quick script to detect paths and define a variable that makes switching between computers easier than commenting-out random lines

detectPaths <- function(){
	# Test if on OSX laptop
	if(dir.exists("~/OneDrive-NIH/") & .Platform$OS.type =="unix"){
		statement <- "I think you're on your laptop"
		comp.type <- "macbookPro"
	}
	else if(dir.exists("/data/CRGGH/heustonef/")){
		statement <- "I think you're on Biowulf/Helix"
		comp.type <- "biowulf"
	}
	else if(dir.exists("/Users/heustonef/OneDrive\ -\ National\ Institutes\ of\ Health")){
		statement <- "I think you're on your work PC and you should probably move to a different system"
		comp.type <- "workPC"
	}
	else{
		print("I'm not sure which computer you're using")
	}
	return(comp.type)
}
