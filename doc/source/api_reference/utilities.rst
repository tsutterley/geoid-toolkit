=========
utilities
=========

Download and management utilities for syncing time and auxiliary files

 - Can list a directory on a ftp host
 - Can list the available static gravity field models from the `ICGEM`__
 - Can download a file from a ftp or http host
 - Checks ``MD5`` or ``sha1`` hashes between local and remote files

.. __: http://icgem.gfz-potsdam.de/home

`Source code`__

.. __: https://github.com/tsutterley/geoid-toolkit/blob/main/geoid_toolkit/utilities.py


General Methods
===============

.. autofunction:: geoid_toolkit.utilities.get_data_path

.. autofunction:: geoid_toolkit.utilities.import_dependency

.. autofunction:: geoid_toolkit.utilities.get_hash

.. autofunction:: geoid_toolkit.utilities.get_git_revision_hash

.. autofunction:: geoid_toolkit.utilities.get_git_status

.. autofunction:: geoid_toolkit.utilities.url_split

.. autofunction:: geoid_toolkit.utilities.convert_arg_line_to_args

.. autofunction:: geoid_toolkit.utilities.get_unix_time

.. autofunction:: geoid_toolkit.utilities.even

.. autofunction:: geoid_toolkit.utilities.ceil

.. autofunction:: geoid_toolkit.utilities.copy

.. autofunction:: geoid_toolkit.utilities.check_ftp_connection

.. autofunction:: geoid_toolkit.utilities.ftp_list

.. autofunction:: geoid_toolkit.utilities.from_ftp

.. autofunction:: geoid_toolkit.utilities._create_default_ssl_context

.. autofunction:: geoid_toolkit.utilities._create_ssl_context_no_verify

.. autofunction:: geoid_toolkit.utilities._set_ssl_context_options

.. autofunction:: geoid_toolkit.utilities.check_connection

.. autofunction:: geoid_toolkit.utilities.from_http

.. autofunction:: geoid_toolkit.utilities.icgem_list
