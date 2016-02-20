HYPHYDIR = "placeholder/batchfiles/";
WDIR="placeholder/";
datafile="temp.fasta";
treefile="temp.tre";
nucfitfile="nuc.fit";
outfile="fel.txt";

inputRedirect = {};
inputRedirect["01"]="Universal";        // Genetic Code
inputRedirect["02"]="New Analysis";     // New analysis
inputRedirect["03"]=WDIR+datafile;      // Alignment file
inputRedirect["04"]="Custom";          // Use custom nucleotide model
inputRedirect["05"]="012345";          // GTR nucleotide model
inputRedirect["06"]=WDIR+treefile;      // Tree file
inputRedirect["07"]=WDIR+nucfitfile;    // Save nucleotide fit to..
inputRedirect["08"]="Estimate dN/dS only";  // Estimate dN/dS w/out branch corrections
inputRedirect["09"]="One rate FEL";     // 1-rate FEL (dS constant across sites)
inputRedirect["10"]="0.01";             // FPR for LRT
inputRedirect["11"]=WDIR+outfile;       // output file


ExecuteAFile (HYPHYDIR + "QuickSelectionDetection_sjs.bf", inputRedirect);
