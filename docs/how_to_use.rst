====================
How to use pyJASPAR
====================

Once you have installed pyJASPAR, you can load the module and connect to the latest release of JASPAR:

.. code-block:: pycon

    >>> from pyjaspar import jaspardb
    
    #get the JASPAR2020 release object, by default it will get the latest release.    
    >>> jdb_obj = jaspardb(release='JASPAR2020')

    #Fetch motif by ID
    >>> motif = jdb_obj.fetch_motif_by_id('MA0095')
    >>> print(motif)
    TF name YY1
    Matrix ID   MA0095.2
    Collection  CORE
    TF class    C2H2 zinc finger factors
    TF family   More than 3 adjacent zinc finger factors
    Species 9606
    Taxonomic group vertebrates
    Accession   ['P25490']
    Data type used  ChIP-seq
    Medline 18950698
    PAZAR ID    TF0000069
    Matrix:
            0      1      2      3      4      5      6      7      8      9     10     11
    A: 1126.00 6975.00 6741.00 2506.00 7171.00   0.00  11.00  13.00 812.00 867.00 899.00 1332.00
    C: 4583.00   0.00  99.00 1117.00   0.00  12.00   0.00   0.00 5637.00 1681.00 875.00 4568.00
    G: 801.00 181.00 268.00 3282.00   0.00   0.00 7160.00 7158.00  38.00 2765.00 4655.00 391.00
    T: 661.00  15.00  63.00 266.00   0.00 7159.00   0.00   0.00 684.00 1858.00 742.00 880.00

    >>> motifs = jdb_obj.fetch_motifs(
        collection = 'CORE',
        tax_group = ['vertebrates', 'insects'],
        tf_class = 'Homeo domain factors',
        tf_family = ['TALE-type homeo domain factors', 'POU domain factors']
        )
    >>> for motif in motifs:
            pass # do something with the motif

