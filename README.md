# X-EISD: Extended Experimental Inferential Structure Determination
X-EISD is a Bayeian approach to perform Experimental Inferential Structure Determination of ensembles for intrinsically disordered proteins.

--  update: June, 2022


## Installation:
You can install eisd from the source:

    git clone https://github.com/THGLab/X-EISD.git
    cd X-EISD
    pip install -e .

## Dependencies:
The dependencies for X-EISD are only numpy and pandas libraries. A proper installation of the package will also install
the requirements.

## Citation:
Please cite the use of X-EISD as:

    (1) James Lincoff, Mojtaba Haghighatlari, Mickael Krzeminski, Joao Teixeira, Gregory Neal-Gomes, Claudu Gardinaru, Julie Forman-Kay, Teresa Head-Gordon, https://www.nature.com/articles/s42004-020-0323-0
    (2) David H. Brookes, and Teresa Head-Gordon, Experimental Inferential Structure Determination of Ensembles for Intrinsically Disordered Proteins, JACS 138, 2016, 4530-4538 DOI: 10.1021/jacs.6b00351

## Getting Started 
You can follow along the sample_script.py in the repository: 

    python sample_script.py     # first modify sample_script based upon your request 
    
Example experimental data files can be found in the data folder. X-EISD program does not support back calculations, though we show some back calculation examples in the EISD_back_calc notebook. Please refer to the above 2 papers for details on experimental data, back calculation methods and determining uncertainties.
      




