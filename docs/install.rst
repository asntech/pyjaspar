==============
How to Install
==============
pyJASPAR is available on `PyPi <https://pypi.python.org/pypi/pyjaspar>`_, through `Bioconda <https://bioconda.github.io/recipes/pyjaspar/README.html>`_, and source code available on `GitHub <https://github.com/asntech/pyjaspar>`_. If you already have a working installation of Python, the easiest way to install the required Python modules is by installing pyJASPAR using ``pip``. 

If you're setting up Python for the first time, we recommend to install it using the `Conda or Miniconda Python distribution <https://conda.io/docs/user-guide/install/index.html>`_. This comes with several helpful scientific and data processing libraries, and available for platforms including Windows, Mac OSX and Linux.

You can use one of the following ways to install pyJASPAR.


Install uisng Conda
====================
We highly recommend to install pyJASPAR using Conda, this will take care of the dependencies. If you already have Conda or Miniconda installed, go ahead and use the below command.

.. code-block:: bash

	conda install -c bioconda pyjaspar

.. note:: This will install all the dependencies and you are ready to use **pyJASPAR**.

Install using `pip`
===================
You can install pyJASPAR from PyPi using pip.

.. code-block:: bash

	pip install pyjaspar

.. note:: Make sure you're using python v3.6 or latest. 



Install from source
===================
You can install a development version by using ``git`` from our GitHub repository at https://github.com/asntech/pyjaspar. 

.. code-block:: bash

    git clone https://github.com/asntech/pyjaspar.git
    cd pyjaspar
    python setup.py sdist install
