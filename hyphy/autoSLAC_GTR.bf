/* SJS 3/14/15. Script to automate SLAC analyses. Uses the MG94xHKY85.*/

HYPHYDIR = "placeholder/batchfiles/";
WDIR="placeholder/";
datafile="temp.fasta";
treefile="temp.tre";
nucfitfile="nuc.fit";
outfile="slac.txt";
pvaluefile="pvalues.txt";

inputRedirect = {};
inputRedirect["01"]="Universal";        // Genetic Code
inputRedirect["02"]="New Analysis";     // New analysis
inputRedirect["03"]=WDIR+datafile;      // Alignment file
inputRedirect["04"]="Custom";          // custom nucleotide model
inputRedirect["05"]="012345";          // GTR
inputRedirect["06"]=WDIR+treefile;      // Tree file
inputRedirect["07"]=WDIR+nucfitfile;    // Save nucleotide fit to..
inputRedirect["08"]="Estimate dN/dS only";  // Estimate dN/dS w/out branch corrections
inputRedirect["09"]="Single Ancestor Counting";     // SLAC 
inputRedirect["10"]="Full tree";             // Analyze the entire tree
inputRedirect["11"]="Averaged";             // Averaged ambiguities
inputRedirect["12"]="Approximate";             // test statistic
inputRedirect["13"]="0.01";             // alpha
inputRedirect["14"]=WDIR+pvaluefile;             // SJS HACKED OPTION FOR SAVING PVALUES!  
inputRedirect["15"]="Export to File";             // save to file
inputRedirect["16"]=WDIR+outfile;             // outfile name
inputRedirect["17"]="Skip";             // skip the rate class estimator


ExecuteAFile (HYPHYDIR + "QuickSelectionDetection_sjs.bf", inputRedirect);
