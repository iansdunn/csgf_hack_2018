# DOE-CSGF

What are the different versions of the files : 
1) gppKerSeq.cpp - Unoptimized sequential version
2) gppOpenMP3.cpp - OpenMP3.5 version
3) gppMPIOpenMP3.cpp - OpenMP3.5 version
4) gppComplex.cpp - The most optimized version which uses a Complex class for its computations.

Makefile - Modify based on version one is working with.
testSmallJob.pbs - a small test job for titan.  
testMedJob.pbs - a medium test job for titan.  
testBigJob.pbs - a Big test job for titan.  

Submit Job command - qsub job.pbs
    Returns a JobId

Check Job Status - showq | grep JobId
