============
utilities.py
============

Download and management utilities for syncing time and auxiliary files

 - Can list a directory on a ftp host
 - Can list the available static gravity field models from the `ICGEM`__
 - Can download a file from a ftp or http host
 - Checks MD5 hashes between local and remote files

.. __: http://icgem.gfz-potsdam.de/home

`Source code`__

.. __: https://github.com/tsutterley/geoid-toolkit/blob/main/geoid_toolkit/utilities.py


General Methods
===============


.. method:: geoid_toolkit.utilities.get_data_path(relpath)

    Get the absolute path within a package from a relative path

    Arguments:

        `relpath`: local relative path as list or string


.. method:: geoid_toolkit.utilities.get_hash(local)

    Get the MD5 hash value from a local file or BytesIO object

    Arguments:

        `local`: BytesIO object or path to file


.. method:: geoid_toolkit.utilities.url_split(s)

    Recursively split a url path into a list

    Arguments:

        `s`: url string


.. method:: geoid_toolkit.utilities.get_unix_time(time_string, format='%Y-%m-%d %H:%M:%S')

    Get the Unix timestamp value for a formatted date string

    Arguments:

        `time_string`: formatted time string to parse

    Keyword arguments:

        `format`: format for input time string


.. method:: geoid_toolkit.utilities.copy(source, destination, verbose=False, move=False)

    Copy or move a file with all system information

    Arguments:

        `source`: source file

        `destination`: copied destination file

    Keyword arguments:

        `verbose`: print file transfer information

        `move`: remove the source file


.. method:: geoid_toolkit.utilities.ftp_list(HOST,timeout=None,basename=False,pattern=None,sort=False)

    List a directory on a ftp host

    Arguments:

        `HOST`: remote ftp host path split as list

    Keyword arguments:

        `timeout`: timeout in seconds for blocking operations

        `basename`: return the file or directory basename instead of the full path

        `pattern`: regular expression pattern for reducing list

        `sort`: sort output list

    Returns:

        `output`: list of items in a directory

        `mtimes`: list of last modification times for items in the directory


.. method:: geoid_toolkit.utilities.from_ftp(HOST,timeout=None,local=None,hash='',chunk=8192,verbose=False,fid=sys.stdout,mode=0o775)

    Download a file from a ftp host

    Arguments:

        `HOST`: remote ftp host path split as list

    Keyword arguments:

        `timeout`: timeout in seconds for blocking operations

        `local`: path to local file

        `hash`: MD5 hash of local file

        `chunk`: chunk size for transfer encoding

        `verbose`: print file transfer information

        `fid`: open file object to print if verbose

        `mode`: permissions mode of output local file


.. method:: geoid_toolkit.utilities.check_connection(HOST)

    Check internet connection

    Arguments:

        `HOST`: remote http host


.. method:: geoid_toolkit.utilities.from_http(HOST,timeout=None,context=ssl.SSLContext(),local=None,hash='',chunk=16384,verbose=False,fid=sys.stdout,mode=0o775)

    Download a file from a http host

    Arguments:

        `HOST`: remote http host path split as list

    Keyword arguments:

        `timeout`: timeout in seconds for blocking operations

        `context`: SSL context for url opener object

        `local`: path to local file

        `hash`: MD5 hash of local file

        `chunk`: chunk size for transfer encoding

        `verbose`: print file transfer information

        `fid`: open file object to print if verbose

        `mode`: permissions mode of output local file


.. method:: geoid_toolkit.utilities.icgem_list(host='http://icgem.gfz-potsdam.de/tom_longtime',timeout=None,parser=lxml.etree.HTMLParser())

    Parse table of gravity field models on the `GFZ International Centre for Global Earth Models (ICGEM)`__ server

    Keyword arguments:

        `host`: url for the GFZ ICGEM gravity field table

        `timeout`: timeout in seconds for blocking operations

        `parser`: HTML parser for lxml

    Returns:

        `colfiles`: dictionary of static file urls mapped by field name

    .. __: http://icgem.gfz-potsdam.de/
