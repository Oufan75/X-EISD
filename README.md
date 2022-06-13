# X-EISD: Extended Experimental Inferential Structure Determination
X-EISD is a Bayeian approach to perform Experimental Inferential Structure Determination of ensembles for intrinsically disordered proteins.

--  update: June 13, 2022


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
You can either follow along the sample_script.py in the repository or use the commnd line interface to run eisd: 

    python sample_script.py     # first modify sample_script based upon your request 
      
or

    eisdshell -d/--datapath, -m/--mode, -s/--structure, -e/--epochs, -o/--output



