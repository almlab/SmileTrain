### Find the database
You'll need to find the file that links the Greengenes IDs with the taxonomies. For example, `gg_13_5/gg_13_5_taxonomy.txt`.

Searching the plain text file is slow, so you can create a fast look-up hash table using the Smile Train script `util/setup_tools/pickle_taxonomy.py`. I recommend placing the pickled file in, e.g., `gg_13_5/pickles/gg_13_5_taxonomy.pkl`.

### Getting the taxonomies
You can use the Smile Train script `util/get_taxonomies.py`. For example, to get the Greengenes taxonomies for an OTU table `otu.txt` that has Greengenes IDs and put them in `tax.txt`, you would call

`/path/to/SmileTrain/util/get_taxonomies.py otu.txt /path/to/gg_13_5/pickles/gg_13_5_taxonomy.pkl > tax.txt`

If you are using a pickled dictionary, this is fast, so you can just run it on the head node of a compute cluster. If you have not pickled the taxonomy list, it would be impolite to run on the head node.