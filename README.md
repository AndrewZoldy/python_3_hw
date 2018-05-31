# python_3_hw


Given script show you spectres of k-mers in your .fastq file, and prints size of genome included in current fastq.

## Launch arguments

Script requires path to input file, you have to define it with **-f** (**--file**) argument:

*-f /path/to/input/fastq_file.fastq*

*--file /path/to/input/fastq_file.fastq*

Also you can set such parametres as : k-mer size **-s** (**--size**), quality threshold **-q** (**--quality**) and noise threshold **-n** (**--noise**) (noise is count of frequences of k-mers under given threshold). 

For this parametres simply use:

*-s 30*

*--size 30*

------

*-q 10*

*--quality 10*

------

*-n 5*

*--noise 5*

By default:
  * K-mer size = 15
  * Quality = 20
  * Noise = 5
  
  
## Results
As result you will obtain 2 plots (without drawen threshold and with it) and value which indicating genome size.
