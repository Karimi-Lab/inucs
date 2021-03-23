# iNucs

Given a nucleosomes file and a DNA-DNA interactions files produced by the [Pairx](https://github.com/4dn-dcic/pairix) program, the `inucs` command line tool bins interactions falling into different nucleosomes and counts them.

[TOC]


------
## Installation

1. Install [Anaconda](https://www.anaconda.com/products/individual) with python 3.8 or above
2. Install the dependencies using Anaconda. In a terminal where anaconda is activated type:

```bash
conda install numpy pandas bokeh
```
3. Then the fully self-contained python script, [`inucs.py`](./inucs.py), may be directly executed; for example: 

```bash
./inucs.py --help
```




------
## Usage



<img src="docs/inucs_workflow.jpg" alt="iNucs Pipeline" style="zoom:40%;" />



As depicted in the figure above, there are two main steps for using `inucs`, namely: `prepare` and `plot`. There is built-in help provided for each stage of the program, which are accessible via command line flags `-h` or `--help`.

Overall help:

```bash
./inucs.py --help
```

which outputs:

```
usage: inucs.py [-h] [-q] {prepare,plot} ...

(omitted for brevity...)
```

The curly brackets `{}` denote that either of `prepare` or `plot` can be used as commands to the program. The `-h` flag is to print this help message, and `-q` suppresses the program progress output.

### The `prepare` Command

In this step, a given potentially large DNA interaction file is broken into smaller pieces and corresponding nucleosome-nucleosome interaction matrices are build.

To access the built-in help for `prepare` command, issue:

```bash
./inucs.py prepare --help
```

which outputs:

```
usage: inucs.py prepare [-h] [-d <working_dir>] [--refresh] [-z] <chroms> <nucs> <interacts>

(omitted for brevity...)
```



  * Input

    All input files can optionally be compressed using the `gzip` format, in which case they need to use `.gz` file extension.

    * `<chroms>`: a file listing the chromosomes of interest

    * `<nucs>`: a file containing nucleosomes

    * `<interacts>`: a file containing DNA-DNA interactions produced by the [Pairx](https://github.com/4dn-dcic/pairix) program

      

  * Output
    
    * `<working_dir>`: a folder containing all the intermediary and cached files. These include the resulting matrices with nucleosome-nucleosome interaction counts. Specifying the `<working_dir>` is optional, and if missing, the program will auto-generate a name based on the name of the interactions input file, `<interacts>`.
    
      If the optional `-z` flag is used, the intermediary files in the `<working_dir>` will be compressed using the `gzip` format.
    
      If the optional `--refresh` flag is used, it starts from scratch and recreates the intermediary files as needed.



### The `plot` Command

The results produced by the `prepare` command are stored in the `<working_dir>` folder, from which the  `plot` command can produce heatmap plots of nucleosome-nucleosome interaction counts matrix. The user specifies the DNA region of interest, and the program finds out nucleosomes within that region, and selects the the submatrix nucleosome interactions existing within the user-specified DNA region. The final output of the `plot` command is a [Bokeh](https://bokeh.org/) interactive plot in `html` format, which can be opened in any standard browser such as Chrome.

To access the built-in help for `plot` command, use:
```bash
./inucs.py plot --help
```

which outputs:
```
usage: inucs.py plot [-h] [-p <outfile_prefix>] [-s] <working_dir> <chrom> <start_region> <end_region>

(omitted for brevity...)
```



* Input

  * `<working_dir>` which is created using the `prepare` command (see above)
  * `<chrom>` such as `III` or `chromX`
  * `<start_region>`  and `<end_region> ` which specify the beginning and end of the region of interest, such as `50000 60000`, within`<chrom>`

* Output

  * The resulting heatmap plot file is also put inside the `<working_dir>` and starts with prefix `plot` or  `<outfile_prefix>` if that is specified by the user. The following pattern is used to generate the resulting plot file name:

    `<working_dir>/<outfile_prefix>_<chrom>_<start_region>_<end_region>.html`

    For example:

    `yeast_wd/plot_I_50000_60000.html`



### Examples

Let us see some example to help make this more clear. As mentioned above `./inucs.py prepare --help` gives the following usage message:

```
usage: inucs.py prepare [-h] [-d <working_dir>] [--refresh] [-z] <chroms> <nucs> <interacts>
```

Thus, as the first step, we can preprocess the input data using the following command:

```bash
./inucs.py prepare chromosomes.tsv nucleosomes.tsv yeast.tsv -d yeast_wd
```

Depending on the size of input files and the system specifications, the `prepare` command may take a few minutes or many hours or days to complete.

The next step is to produce plots, which we saw its usage above, again using the built-in help `./inucs.py plot --help`:

```
usage: inucs.py plot [-h] [-p <outfile_prefix>] [-s] <working_dir> <chrom> <start_region> <end_region>
```

Thus, continuing after the `prepare` command above, to get interactions on chromosome `I` in the region between 50000 and 60000, we can use:


```bash
./inucs.py plot yeast_wd I 50000 60000
```



<img src="docs/plot_I_50000_60000.jpg" alt="plot_I_50000_60000" style="zoom:35%;" />

If the figure does not show automatically you can open it manually. 
For example, on a **macOS** system, you can use Find to navigate to `yeast_wd`  working directory and them open the plot file. Alternatively, you can simply use appropriate commands in terminal to open the plot file. For example, in **macOS** use:

```bash
open yeast_wd/plot_I_50000_60000.html
```

Or, on **Windows**, you may use:

```bash
explorer.exe yeast_wd\plot_I_50000_60000.html
```


Note that if you are using a **Linux** terminal in **Windows** using [WSL](https://docs.microsoft.com/en-us/windows/wsl/install-win10), then you may need to escape the backslash characters, i.e., use `\\` instead of `\`.



------
## Implementation

### Efficiency

Given that the input data for `inucs` is expected to be very large, it is important to take time complexity of the underlying algorithms very seriously. Time complexity models the expected amount of time needed for an algorithm to compete a task in terms of the size of its input. For example, a sorting algorithm may take as input *n* numbers and it may take a time proportional to *n log n* to sort the array. We denote such performance or time complexity using the standard big-*O* notation as follows:
$$
O(n\ log\ n)
$$
There has been considerable work put into sorting algorithms, and the above time complexity is the state of the art for many important search algorithms.

For `inucs`, we have managed to reduce the problem of matching DNA interactions with nucleosomes into a sorting problem. This in turns allows us to leverage the power of existing sorting algorithms to optimize to solution to our problem of nucleosome interactions.

Another benefit of reducing the nucleosome problem into existing frameworks is that now `inucs` utilizes of vectorization capabilities provided by NumPy/Pandas with hardware support. This is all without adding extra visible complications within `inucs` to add vectorization, or in other words, without reinventing the wheel.

We are planning for the next version of `inucs` to add support for parallelism for higher utilization of modern multicore CPUs.


### Scalability
We have put a large effort to make `inucs` scalable as much as possible given our limited development time. In its correct state, `inucs` can handle large amounts of input data without requiring exceeding computational resources. For example, `inucs` can process more than 3.2 billion human chromosome interactions (over 261 GB file size), while running a regular PC as long as it has 40 GB of RAM or more. (It took about 16 hours to complete in our testing.) 

Scalability in `inucs` is achieved primarily by breaking down the data into smaller pieces in different stages of running the program. Therefore, at any given time, there is only some manageable chunk of data in memory, and the intermediary results are constantly written on storage space and read back in as needed.

Some important breakdowns of the data are as follows
* Input: Interactions between chromosome locations

   Broken into chromosomes and orientations (strands ++, --, +-, -+); e.g., for human data that would be 24 chromosomes times 4 orientations, which means 96 data files extracted from original input.
   
* Nucleosome-Nucleosome Interaction Matrix

   Broken into chromosomes and orientations (e.g., 96 files similar to above)
   
* Submatrix selection for plotting

   A smaller range from the above matrix is selected by specifying: chromosome beginning-of-range end-of-range


### Runtime Measurements

As mentioned, currently we have not taken advantage of multiprocessing yet for the program. Thus, our runtime measurements are using a single CPU core on both the PC and HPC server examples below.

#### Example 1: Yeast Data

The `inucs` application can process yeast data on a regular PC. For one example run, we have used a

> Laptop, using one core from a Intel(R) Core(TM) i7-8565U CPU, with 16 GB of RAM, running Windows 10.

The yeast data used:

* For `<nucs>`, 77,060 nucleosomes
* For `<interacts>`, 24,163,427 DNA interactions
* With a resulting nucleosome interaction **matrix size** of **77,060**nucs x **77,060**nucs 

The `prepare` commands, takes less than ***4 minutes*** to complete for this example. 

The time for `plot` is less than *15 seconds*.

#### Example 2: Human Data

The application takes much longer to process human data as expected. For this example run, we have used an

> HPC server, using one core from an AMD EPYC 7552 CPU, with 40 GB of RAM, running CentOS Linux 7.6.

The human data used:

* For `<nucs>`, 13,811,032 nucleosomes
* For `<interacts>`, 3,220,503,431 DNA interactions
* With a resulting nucleosome interaction **matrix size** of **13,811,032**nucs x **13,811,032**nucs

The `prepare` commands, takes about ***980 minutes*** (or less than 17 hours) to complete for this example. One important note here is that most of the time is spent to break down the large input file into smaller manageable files, each of which contain data only for one chromosome and one orientation (i.e., --, ++, -+, or +-). More, specifically it takes about 700 minutes (less that 12 hours) to break down the file, and the remaining 280 minutes (less that 5  hours) to complete the interaction matrix calculations.

The time for `plot` is about than *3.5 minutes*.
