This ``bfdb`` branch:
  * contains the longstanding history of the repo, excluding that developed in Psi4
  * is relevant to ``The BioFragment Database (BFDb): An Open-Data Platform for Computational Chemistry Analysis of Noncovalent Interactions,'' L. A. Burns, J. C. Faver, Z. Zheng, M. S. Marshall, D. G. A. Smith, K. Vanommeslaeghe, A. D. MacKerell, K. M. Merz, and C. D. Sherrill, J. Chem. Phys. 147, 161727 (2017) (http://aip.scitation.org/doi/10.1063/1.5001028)
    

.. image:: https://zenodo.org/badge/5744708.svg
   :target: https://zenodo.org/badge/latestdoi/5744708
   
qcdb
====

Quantum Chemistry Common Driver and Databases
---------------------------------------------

Modules and scripts for managing databases of chemical systems and storing
quantum chemical information. Promoting pandas for storage and
manipulation of quantum chemical information. Facilitating
interoperability among quantum chemistry codes. Closely tied in with the
Psi4 project: www.psicode.org and https://github.com/psi4 .

Because one large file is under git lfs control, you may need to issue the following after cloning. Git itself must be >=1.8.5::

    conda install git-lfs -c conda-forge
    git lfs fetch
    git lfs checkout
