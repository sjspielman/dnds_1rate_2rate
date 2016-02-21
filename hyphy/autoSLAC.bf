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
inputRedirect["04"]="Default";          // Use HKY85 and MG94xHKY85.
inputRedirect["05"]=WDIR+treefile;      // Tree file
inputRedirect["06"]=WDIR+nucfitfile;    // Save nucleotide fit to..
inputRedirect["07"]="Estimate dN/dS only";  // Estimate dN/dS w/out branch corrections
inputRedirect["08"]="Single Ancestor Counting";     // SLAC 
inputRedirect["09"]="Full tree";             // Analyze the entire tree
inputRedirect["10"]="Averaged";             // Averaged ambiguities
inputRedirect["11"]="Approximate";             // test statistic
inputRedirect["12"]="0.01";             // alpha
inputRedirect["13"]=WDIR+pvaluefile;             // SJS HACKED OPTION FOR SAVING PVALUES!  
inputRedirect["14"]="Export to File";             // save to file
inputRedirect["15"]=WDIR+outfile;             // outfile name
inputRedirect["16"]="Skip";             // skip the rate class estimator


ExecuteAFile (HYPHYDIR + "QuickSelectionDetection_sjs.bf", inputRedirect);
