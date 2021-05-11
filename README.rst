pyJASPAR
--------

    A Pythonic interface to JASPAR transcription factor motifs

**pyJASPAR** uses *Biopython* and *SQLite3* to provide a serverless interface to `JASPAR database <http://jaspar.genereg.net>`_ to query and access TF motif profiles across various releases of JASPAR.


.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.4509415.svg
   :target: https://doi.org/10.5281/zenodo.4509415

.. image:: https://travis-ci.org/asntech/pyjaspar.svg?branch=main
    :target: https://travis-ci.org/asntech/pyjaspar

.. image:: https://img.shields.io/pypi/pyversions/pyjaspar.svg
    :target: https://www.python.org

.. image:: https://img.shields.io/pypi/v/pyjaspar.svg
    :target: https://pypi.python.org/pypi/pyjaspar

.. image:: https://anaconda.org/bioconda/pyjaspar/badges/version.svg
    :target: https://anaconda.org/bioconda/pyjaspar

.. image:: https://anaconda.org/bioconda/pyjaspar/badges/downloads.svg
    :target: https://bioconda.github.io/recipes/pyjaspar/README.html

.. image:: https://img.shields.io/github/issues/asntech/pyjaspar.svg
    :target: https://github.com/asntech/pyjaspar/issues


pyJASPAR provides access to the following releases of JASPAR database:
    - `JASPAR2020 <http://jaspar2020.genereg.net>`_
    - `JASPAR2018 <http://jaspar2018.genereg.net>`_
    - `JASPAR2016 <http://jaspar2016.genereg.net>`_
    - `JASPAR2014 <http://jaspar2014.genereg.net>`_


**Note**: This is a serverless SQLite wrapper around the Biopython JASPAR module `Bio.motifs.jaspar.db` which requires JASPAR MySQL database sever connection details.


Documentation
-------------

**A detailed documentation is available in different formats:**  `HTML <http://pyjaspar.readthedocs.org>`_ | `PDF <http://readthedocs.org/projects/pyjaspar/downloads/pdf/latest/>`_ | `ePUB <http://readthedocs.org/projects/pyjaspar/downloads/epub/latest/>`_


Installation
------------

Quick installation using conda
================================
pyJASPAR is available on `Bioconda <https://anaconda.org/bioconda/pyjaspar>`_ for installation via ``conda``.

.. code-block:: bash

	conda install -c bioconda pyjaspar


Install using pip
==================
pyJASPAR is also available on `PyPi <https://pypi.org/project/pyjaspar/>`_ for installation via ``pip``.

.. code-block:: bash

	pip install pyjaspar
	

pyJASPAR uses BioPython and it supports python ``3.x``. 

Install pyjaspar from source
=============================
You can install a development version by using ``git`` from GitHub.


Install development version from `GitHub`
==========================================
If you have `git` installed, use this:

.. code-block:: bash

    git clone https://github.com/asntech/pyjaspar.git
    cd pyjaspar
    python setup.py sdist install

How to use pyJASPAR?
--------------------

Once you have installed pyjaspar, you can create jaspardb class object:

.. code-block:: pycon

    >>> from pyjaspar import jaspardb
    
    #Connect the JASPAR2020 release object    
    >>> jdb_obj = jaspardb(release='JASPAR2020')

    #Fetch motif by ID
    >>> motif = jdb_obj.fetch_motif_by_id('MA0095.2')
    >>> print(motif.name)
    YY1

    #Fetch motifs by TF name
    >>> motifs = jdb_obj.fetch_motifs_by_name('CTCF')
    >>> print(len(motifs))
    2

    #Get CORE vertebrates collection
    >>> motifs = jdb_obj.fetch_motifs(
        collection = 'CORE',
        tax_group = ['vertebrates']
        )
    >>> print(len(motifs))
    746
    ## loop through the motifs list and perform analysis
    >>> for motif in motifs:
            pass

**Note**: Above methods return `Bio.motifs.jaspar.Motif` object. You can find more details `here <http://biopython.org/DIST/docs/tutorial/Tutorial.html#sec262>`_ 


Find available releases
=======================
.. code-block:: pycon
    
    >>> print(jdb_obj.get_releases())
    Available JASPAR releases are: ['JASPAR2020', 'JASPAR2018', 'JASPAR2016', 'JASPAR2014']


Cite
=====
    - Aziz Khan. pyJASPAR: a Pythonic interface to JASPAR transcription factor motifs. (2021). doi:10.5281/zenodo.4509415

.. code-block:: bash

    @software{aziz_khan_2021_4509415,
      author       = {Aziz Khan},
      title        = {{pyJASPAR: a Pythonic interface to JASPAR 
                       transcription factor motifs}},
      month        = feb,
      year         = 2021,
      publisher    = {Zenodo},
      version      = {v1.5.5},
      doi          = {10.5281/zenodo.4509415},
      url          = {https://doi.org/10.5281/zenodo.4509415}
    }
