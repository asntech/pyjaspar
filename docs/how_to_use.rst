====================
How to use pyJASPAR
====================

Once you have installed pyJASPAR, you can load the module and connect to the latest release of JASPAR:

.. code-block:: pycon

    >>> from pyjaspar import jaspardb

Connect to the version of JASPAR you're interested in. This will return jaspardb class object.
For example here we're getting the JASPAR2018.

.. code-block:: pycon

    >>> jdb_obj = jaspardb(release='JASPAR2018')

You can also check JASPAR version you are connected to using:

.. code-block:: pycon

    >>> print(jdb_obj.release)
    JASPAR2018

By default it is set to latest release/version of JASPAR database. For example.

.. code-block:: pycon

    >>> jdb_obj = jaspardb()
    >>> print(jdb_obj.release)
    JASPAR2020


Get available releases
----------------------
You can find the available releases/version of JASPAR using.


.. code-block:: pycon

    >>> print(jdb_obj.get_releases())
    ['JASPAR2020', 'JASPAR2018', 'JASPAR2016', 'JASPAR2014']


Get motif by using JASPAR ID
----------------------------
If you want to get the motif details for a specific TF using the JASPAR ID. If you skip the version of motif, it will return the latest version. 

.. code-block:: pycon

    >>> motif = jdb_obj.fetch_motif_by_id('MA0095.2')

Printing the motif will all the associated meta-information stored in the JASPAR database cluding the matric counts.


.. code-block:: pycon

    >>> print(motif)
    TF name	YY1
	Matrix ID	MA0095.2
	Collection	CORE
	TF class	C2H2 zinc finger factors
	TF family	More than 3 adjacent zinc finger factors
	Species	9606
	Taxonomic group	vertebrates
	Accession	['P25490']
	Data type used	ChIP-seq
	Medline	18950698
	Matrix:
	        0      1      2      3      4      5      6      7      8      9     10     11
	A: 1126.00 6975.00 6741.00 2506.00 7171.00   0.00  11.00  13.00 812.00 867.00 899.00 1332.00
	C: 4583.00   0.00  99.00 1117.00   0.00  12.00   0.00   0.00 5637.00 1681.00 875.00 4568.00
	G: 801.00 181.00 268.00 3282.00   0.00   0.00 7160.00 7158.00  38.00 2765.00 4655.00 391.00
	T: 661.00  15.00  63.00 266.00   0.00 7159.00   0.00   0.00 684.00 1858.00 742.00 880.00


Get the count matrix using `.counts`


.. code-block:: pycon

    >>> print(motif.counts)
            0      1      2      3      4      5      6      7      8      9     10     11
	A: 1126.00 6975.00 6741.00 2506.00 7171.00   0.00  11.00  13.00 812.00 867.00 899.00 1332.00
	C: 4583.00   0.00  99.00 1117.00   0.00  12.00   0.00   0.00 5637.00 1681.00 875.00 4568.00
	G: 801.00 181.00 268.00 3282.00   0.00   0.00 7160.00 7158.00  38.00 2765.00 4655.00 391.00
	T: 661.00  15.00  63.00 266.00   0.00 7159.00   0.00   0.00 684.00 1858.00 742.00 880.00


Get motifs by TF name
-----------------------
You can use the `fetch_motifs_by_name` function to find motifs by TF name. This method returns a list of motifs for the same TF name across taxonomic group. For example, below search will return two CTCF motifs one in vertebrates and another in plants taxon.

.. code-block:: pycon

    >>> motifs = jdb_obj.fetch_motifs_by_name("CTCF")
    >>> print(len(motifs))
    2
    >>> print(motifs)
    TF name	CTCF
	Matrix ID	MA0139.1
	Collection	CORE
	TF class	C2H2 zinc finger factors
	TF family	More than 3 adjacent zinc finger factors
	Species	9606
	Taxonomic group	vertebrates
	Accession	['P49711']
	Data type used	ChIP-seq
	Medline	17512414
	Matrix:
	        0      1      2      3      4      5      6      7      8      9     10     11     12     13     14     15     16     17     18
	A:  87.00 167.00 281.00  56.00   8.00 744.00  40.00 107.00 851.00   5.00 333.00  54.00  12.00  56.00 104.00 372.00  82.00 117.00 402.00
	C: 291.00 145.00  49.00 800.00 903.00  13.00 528.00 433.00  11.00   0.00   3.00  12.00   0.00   8.00 733.00  13.00 482.00 322.00 181.00
	G:  76.00 414.00 449.00  21.00   0.00  65.00 334.00  48.00  32.00 903.00 566.00 504.00 890.00 775.00   5.00 507.00 307.00  73.00 266.00
	T: 459.00 187.00 134.00  36.00   2.00  91.00  11.00 324.00  18.00   3.00   9.00 341.00   8.00  71.00  67.00  17.00  37.00 396.00  59.00


	TF name	CTCF
	Matrix ID	MA0531.1
	Collection	CORE
	TF class	C2H2 zinc finger factors
	TF family	More than 3 adjacent zinc finger factors
	Species	7227
	Taxonomic group	insects
	Accession	['Q9VS55']
	Data type used	ChIP-chip
	Medline	17616980
	Matrix:
	        0      1      2      3      4      5      6      7      8      9     10     11     12     13     14
	A: 306.00 313.00 457.00 676.00 257.00 1534.00 202.00 987.00   2.00   0.00   2.00 124.00   1.00  79.00 231.00
	C: 876.00 1147.00 383.00 784.00 714.00   1.00   0.00   0.00   4.00   0.00   0.00 1645.00   0.00 1514.00 773.00
	G: 403.00 219.00 826.00 350.00  87.00 192.00 1700.00 912.00 311.00 1902.00 1652.00   3.00 1807.00   8.00 144.00
	T: 317.00 223.00 236.00  92.00 844.00 175.00   0.00   3.00 1585.00   0.00 248.00 130.00  94.00 301.00 754.00


Search motifs based on meta-info
---------------------------------
A more commonly used function is `fetch_motifs` helps you to get motifs which match a specified set of criteria.
You can query the database based on the available meta-information in the database.

For example, here we are gettting the widely used CORE collection for vertebrates. It returns a list of non-redundent motifs. 

.. code-block:: pycon

    >>> motifs = jdb_obj.fetch_motifs(
    collection = 'CORE',
    tax_group = ['vertebrates']
    )
    >>> print(len(motifs))
    746

You can loop through these motifs and perform your analysis.

.. code-block:: pycon

    >>> for motif in motifs:
    		print(motif.matrix_id)
    MA0004.1
	MA0006.1
	-
	-
	-
	MA0528.2
	MA0609.2