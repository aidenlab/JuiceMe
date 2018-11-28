# JuiceMe
Hi-Culfite Processing Pipeline

JuiceMe is heavily based on the <a href="https://github.com/aidenlab/juicer">Juicer</a> pipeline.
Please review the <a href="https://github.com/aidenlab/juicer/wiki">Juicer documentation</a> first
to understand how the JuiceMe pipeline works.

# Dependencies
In addition to the <a href="https://github.com/aidenlab/juicer/wiki/Installation#dependencies">Juicer dependencies</a>,
you will need:

* <a href="https://github.com/brentp/bwa-meth">bwa-meth</a>
* <a href="http://www.htslib.org/">Samtools</a>
* Python
* <a href="https://github.com/dpryan79/MethylDackel">MethylDackel</a>

Be sure to download the latest <a href="https://github.com/aidenlab/juicer/wiki/Download">Juicer Tools jar</a> 
and to put the MethylDackel executable in the scripts directory. 

You will also need to create index files using bwa-meth and put these in your references directory.

Otherwise, the pipeline proceeds as described in the Juicer documentation and the setup is the same.

# Analysis
Python scripts for the various methylation analyses are included in the Analysis folder.


