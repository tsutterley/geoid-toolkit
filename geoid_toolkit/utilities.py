#!/usr/bin/env python
u"""
utilities.py
Written by Tyler Sutterley (05/2023)
Download and management utilities for syncing time and auxiliary files

PYTHON DEPENDENCIES:
    lxml: processing XML and HTML in Python
        https://pypi.python.org/pypi/lxml

UPDATE HISTORY:
    Updated 05/2023: use pathlib to define and operate on paths
    Updated 01/2023: add default ssl context attribute with protocol
    Updated 12/2022: functions for managing and maintaining git repositories
    Updated 04/2022: updated docstrings to numpy documentation format
    Updated 10/2021: using python logging for handling verbose output
    Updated 09/2021: added generic list from Apache http server
    Updated 07/2021: added unique filename opener for log files
    Updated 06/2021: add parser for converting file lines to arguments
    Updated 03/2021: added sha1 option for retrieving file hashes
    Updated 11/2020: added list function for finding files on the GFZ ICGEM
    Updated 09/2020: copy from http and https to bytesIO object in chunks
    Written 08/2020
"""
from __future__ import print_function, division, annotations

import sys
import os
import re
import io
import ssl
import ftplib
import shutil
import socket
import inspect
import hashlib
import logging
import pathlib
import dateutil
import warnings
import posixpath
import lxml.etree
import subprocess
import calendar,time
if sys.version_info[0] == 2:
    import urllib2
else:
    import urllib.request as urllib2

# PURPOSE: get absolute path within a package from a relative path
def get_data_path(relpath: list | str | pathlib.Path):
    """
    Get the absolute path within a package from a relative path

    Parameters
    ----------
    relpath: list, str or pathlib.Path
        relative path
    """
    # current file path
    filename = inspect.getframeinfo(inspect.currentframe()).filename
    filepath = pathlib.Path(filename).absolute().parent
    if isinstance(relpath, list):
        # use *splat operator to extract from list
        return filepath.joinpath(*relpath)
    elif isinstance(relpath, str):
        return filepath.joinpath(relpath)

# PURPOSE: get the hash value of a file
def get_hash(
        local: str | io.IOBase | pathlib.Path,
        algorithm: str = 'MD5'
    ):
    """
    Get the hash value from a local file or ``BytesIO`` object

    Parameters
    ----------
    local: obj, str or pathlib.Path
        BytesIO object or path to file
    algorithm: str, default 'MD5'
        hashing algorithm for checksum validation

            - ``'MD5'``: Message Digest
            - ``'sha1'``: Secure Hash Algorithm
    """
    # check if open file object or if local file exists
    if isinstance(local, io.IOBase):
        if (algorithm == 'MD5'):
            return hashlib.md5(local.getvalue()).hexdigest()
        elif (algorithm == 'sha1'):
            return hashlib.sha1(local.getvalue()).hexdigest()
    elif isinstance(local, (str, pathlib.Path)):
        # generate checksum hash for local file
        local = pathlib.Path(local).expanduser()
        # if file currently doesn't exist, return empty string
        if not local.exists():
            return ''
        # open the local_file in binary read mode
        with local.open(mode='rb') as local_buffer:
            # generate checksum hash for a given type
            if (algorithm == 'MD5'):
                return hashlib.md5(local_buffer.read()).hexdigest()
            elif (algorithm == 'sha1'):
                return hashlib.sha1(local_buffer.read()).hexdigest()
    else:
        return ''

# PURPOSE: get the git hash value
def get_git_revision_hash(
        refname: str = 'HEAD',
        short: bool = False
    ):
    """
    Get the ``git`` hash value for a particular reference

    Parameters
    ----------
    refname: str, default HEAD
        Symbolic reference name
    short: bool, default False
        Return the shorted hash value
    """
    # get path to .git directory from current file path
    filename = inspect.getframeinfo(inspect.currentframe()).filename
    basepath = pathlib.Path(filename).absolute().parent.parent
    gitpath = basepath.joinpath('.git')
    # build command
    cmd = ['git', f'--git-dir={gitpath}', 'rev-parse']
    cmd.append('--short') if short else None
    cmd.append(refname)
    # get output
    with warnings.catch_warnings():
        return str(subprocess.check_output(cmd), encoding='utf8').strip()

# PURPOSE: get the current git status
def get_git_status():
    """Get the status of a ``git`` repository as a boolean value
    """
    # get path to .git directory from current file path
    filename = inspect.getframeinfo(inspect.currentframe()).filename
    basepath = pathlib.Path(filename).absolute().parent.parent
    gitpath = basepath.joinpath('.git')
    # build command
    cmd = ['git', f'--git-dir={gitpath}', 'status', '--porcelain']
    with warnings.catch_warnings():
        return bool(subprocess.check_output(cmd))

# PURPOSE: recursively split a url path
def url_split(s: str):
    """
    Recursively split a url path into a list

    Parameters
    ----------
    s: str
        url string
    """
    head, tail = posixpath.split(s)
    if head in ('http:','https:','ftp:','s3:'):
        return s,
    elif head in ('', posixpath.sep):
        return tail,
    return url_split(head) + (tail,)

# PURPOSE: convert file lines to arguments
def convert_arg_line_to_args(arg_line):
    """
    Convert file lines to arguments

    Parameters
    ----------
    arg_line: str
        line string containing a single argument and/or comments
    """
    # remove commented lines and after argument comments
    for arg in re.sub(r'\#(.*?)$',r'',arg_line).split():
        if not arg.strip():
            continue
        yield arg

# PURPOSE: returns the Unix timestamp value for a formatted date string
def get_unix_time(
        time_string: str,
        format: str = '%Y-%m-%d %H:%M:%S'
    ):
    """
    Get the Unix timestamp value for a formatted date string

    Parameters
    ----------
    time_string: str
        formatted time string to parse
    format: str, default '%Y-%m-%d %H:%M:%S'
        format for input time string
    """
    try:
        parsed_time = time.strptime(time_string.rstrip(), format)
    except (TypeError, ValueError):
        pass
    else:
        return calendar.timegm(parsed_time)
    # try parsing with dateutil
    try:
        parsed_time = dateutil.parser.parse(time_string.rstrip())
    except (TypeError, ValueError):
        return None
    else:
        return parsed_time.timestamp()

# PURPOSE: output a time string in isoformat
def isoformat(time_string: str):
    """
    Reformat a date string to ISO formatting

    Parameters
    ----------
    time_string: str
        formatted time string to parse
    """
    # try parsing with dateutil
    try:
        parsed_time = dateutil.parser.parse(time_string.rstrip())
    except (TypeError, ValueError):
        return None
    else:
        return parsed_time.isoformat()

# PURPOSE: rounds a number to an even number less than or equal to original
def even(value: float):
    """
    Rounds a number to an even number less than or equal to original

    Parameters
    ----------
    value: float
        number to be rounded
    """
    return 2*int(value//2)

# PURPOSE: rounds a number upward to its nearest integer
def ceil(value: float):
    """
    Rounds a number upward to its nearest integer

    Parameters
    ----------
    value: float
        number to be rounded upward
    """
    return -int(-value//1)

# PURPOSE: make a copy of a file with all system information
def copy(
        source: str | pathlib.Path,
        destination: str | pathlib.Path,
        move: bool = False,
        **kwargs
    ):
    """
    Copy or move a file with all system information

    Parameters
    ----------
    source: str or pathlib.Path
        source file
    destination: str or pathlib.Path
        copied destination file
    move: bool, default False
        remove the source file
    """
    source = pathlib.Path(source).expanduser().absolute()
    destination = pathlib.Path(destination).expanduser().absolute()
    # log source and destination
    logging.info(f'{str(source)} -->\n\t{str(destination)}')
    shutil.copyfile(source, destination)
    shutil.copystat(source, destination)
    # remove the original file if moving
    if move:
        source.unlink()

# PURPOSE: check ftp connection
def check_ftp_connection(
        HOST: str,
        username: str | None = None,
        password: str | None = None
    ):
    """
    Check internet connection with ftp host

    Parameters
    ----------
    HOST: str
        remote ftp host
    username: str or NoneType
        ftp username
    password: str or NoneType
        ftp password
    """
    # attempt to connect to ftp host
    try:
        f = ftplib.FTP(HOST)
        f.login(username, password)
        f.voidcmd("NOOP")
    except IOError:
        raise RuntimeError('Check internet connection')
    except ftplib.error_perm:
        raise RuntimeError('Check login credentials')
    else:
        return True

# PURPOSE: list a directory on a ftp host
def ftp_list(
        HOST: str | list,
        username: str | None = None,
        password: str | None = None,
        timeout: int | None = None,
        basename: bool = False,
        pattern: str | None = None,
        sort: bool = False
    ):
    """
    List a directory on a ftp host

    Parameters
    ----------
    HOST: str or list
        remote ftp host path split as list
    username: str or NoneType
        ftp username
    password: str or NoneType
        ftp password
    timeout: int or NoneType, default None
        timeout in seconds for blocking operations
    basename: bool, default False
        return the file or directory basename instead of the full path
    pattern: str or NoneType, default None
        regular expression pattern for reducing list
    sort: bool, default False
        sort output list

    Returns
    -------
    output: list
        items in a directory
    mtimes: list
        last modification times for items in the directory
    """
    # verify inputs for remote ftp host
    if isinstance(HOST, str):
        HOST = url_split(HOST)
    # try to connect to ftp host
    try:
        ftp = ftplib.FTP(HOST[0],timeout=timeout)
    except (socket.gaierror,IOError):
        raise RuntimeError(f'Unable to connect to {HOST[0]}')
    else:
        ftp.login(username,password)
        # list remote path
        output = ftp.nlst(posixpath.join(*HOST[1:]))
        # get last modified date of ftp files and convert into unix time
        mtimes = [None]*len(output)
        # iterate over each file in the list and get the modification time
        for i,f in enumerate(output):
            try:
                # try sending modification time command
                mdtm = ftp.sendcmd(f'MDTM {f}')
            except ftplib.error_perm:
                # directories will return with an error
                pass
            else:
                # convert the modification time into unix time
                mtimes[i] = get_unix_time(mdtm[4:], format="%Y%m%d%H%M%S")
        # reduce to basenames
        if basename:
            output = [posixpath.basename(i) for i in output]
        # reduce using regular expression pattern
        if pattern:
            i = [i for i,f in enumerate(output) if re.search(pattern,f)]
            # reduce list of listed items and last modified times
            output = [output[indice] for indice in i]
            mtimes = [mtimes[indice] for indice in i]
        # sort the list
        if sort:
            i = [i for i,j in sorted(enumerate(output), key=lambda i: i[1])]
            # sort list of listed items and last modified times
            output = [output[indice] for indice in i]
            mtimes = [mtimes[indice] for indice in i]
        # close the ftp connection
        ftp.close()
        # return the list of items and last modified times
        return (output, mtimes)

# PURPOSE: download a file from a ftp host
def from_ftp(
        HOST: str | list,
        username: str | None = None,
        password: str | None = None,
        timeout: int | None = None,
        local: str | pathlib.Path | None = None,
        hash: str = '',
        chunk: int = 8192,
        verbose: bool = False,
        fid=sys.stdout,
        mode: oct = 0o775
    ):
    """
    Download a file from a ftp host

    Parameters
    ----------
    HOST: str or list
        remote ftp host path
    username: str or NoneType
        ftp username
    password: str or NoneType
        ftp password
    timeout: int or NoneType, default None
        timeout in seconds for blocking operations
    local: str, pathlib.Path or NoneType, default None
        path to local file
    hash: str, default ''
        MD5 hash of local file
    chunk: int, default 8192
        chunk size for transfer encoding
    verbose: bool, default False
        print file transfer information
    fid: obj, default sys.stdout
        open file object to print if verbose
    mode: oct, default 0o775
        permissions mode of output local file

    Returns
    -------
    remote_buffer: obj
        BytesIO representation of file
    """
    # create logger
    loglevel = logging.INFO if verbose else logging.CRITICAL
    logging.basicConfig(stream=fid, level=loglevel)
    # verify inputs for remote ftp host
    if isinstance(HOST, str):
        HOST = url_split(HOST)
    # try downloading from ftp
    try:
        # try to connect to ftp host
        ftp = ftplib.FTP(HOST[0], timeout=timeout)
    except (socket.gaierror,IOError):
        raise RuntimeError(f'Unable to connect to {HOST[0]}')
    else:
        ftp.login(username,password)
        # remote path
        ftp_remote_path = posixpath.join(*HOST[1:])
        # copy remote file contents to bytesIO object
        remote_buffer = io.BytesIO()
        ftp.retrbinary(f'RETR {ftp_remote_path}',
            remote_buffer.write, blocksize=chunk)
        remote_buffer.seek(0)
        # save file basename with bytesIO object
        remote_buffer.filename = HOST[-1]
        # generate checksum hash for remote file
        remote_hash = hashlib.md5(remote_buffer.getvalue()).hexdigest()
        # get last modified date of remote file and convert into unix time
        mdtm = ftp.sendcmd(f'MDTM {ftp_remote_path}')
        remote_mtime = get_unix_time(mdtm[4:], format="%Y%m%d%H%M%S")
        # compare checksums
        if local and (hash != remote_hash):
            # convert to absolute path
            local = pathlib.Path(local).expanduser().absolute()
            # create directory if non-existent
            local.parent.mkdir(mode=mode, parents=True, exist_ok=True)
            # print file information
            args = (posixpath.join(*HOST), str(local))
            logging.info('{0} -->\n\t{1}'.format(*args))
            # store bytes to file using chunked transfer encoding
            remote_buffer.seek(0)
            with local.open(mode='wb') as f:
                shutil.copyfileobj(remote_buffer, f, chunk)
            # change the permissions mode
            local.chmod(mode)
            # keep remote modification time of file and local access time
            os.utime(local, (local.stat().st_atime, remote_mtime))
        # close the ftp connection
        ftp.close()
        # return the bytesIO object
        remote_buffer.seek(0)
        return remote_buffer

# default ssl context
_default_ssl_context = ssl.SSLContext(ssl.PROTOCOL_TLS)

# PURPOSE: check internet connection
def check_connection(HOST: str, context=_default_ssl_context):
    """
    Check internet connection with http host

    Parameters
    ----------
    HOST: str
        remote http host
    context: obj, default ssl.SSLContext(ssl.PROTOCOL_TLS)
        SSL context for ``urllib`` opener object
    """
    # attempt to connect to http host
    try:
        urllib2.urlopen(HOST, timeout=20, context=context)
    except urllib2.URLError as exc:
        raise RuntimeError('Check internet connection') from exc
    else:
        return True

# PURPOSE: list a directory on an Apache http Server
def http_list(
        HOST: str | list,
        timeout: int | None = None,
        context = _default_ssl_context,
        parser = lxml.etree.HTMLParser(),
        format: str = '%Y-%m-%d %H:%M',
        pattern: str = '',
        sort: bool = False
    ):
    """
    List a directory on an Apache http Server

    Parameters
    ----------
    HOST: str or list
        remote http host path
    timeout: int or NoneType, default None
        timeout in seconds for blocking operations
    context: obj, default ssl.SSLContext(ssl.PROTOCOL_TLS)
        SSL context for ``urllib`` opener object
    parser: obj, default lxml.etree.HTMLParser()
        HTML parser for ``lxml``
    format: str, default '%Y-%m-%d %H:%M'
        format for input time string
    pattern: str, default ''
        regular expression pattern for reducing list
    sort: bool, default False
        sort output list

    Returns
    -------
    colnames: list
        column names in a directory
    collastmod: list
        last modification times for items in the directory
    """
    # verify inputs for remote http host
    if isinstance(HOST, str):
        HOST = url_split(HOST)
    # try listing from http
    try:
        # Create and submit request.
        request = urllib2.Request(posixpath.join(*HOST))
        response = urllib2.urlopen(request, timeout=timeout, context=context)
    except (urllib2.HTTPError, urllib2.URLError):
        raise Exception('List error from {0}'.format(posixpath.join(*HOST)))
    else:
        # read and parse request for files (column names and modified times)
        tree = lxml.etree.parse(response, parser)
        colnames = tree.xpath('//tr/td[not(@*)]//a/@href')
        # get the Unix timestamp value for a modification time
        collastmod = [get_unix_time(i,format=format)
            for i in tree.xpath('//tr/td[@align="right"][1]/text()')]
        # reduce using regular expression pattern
        if pattern:
            i = [i for i,f in enumerate(colnames) if re.search(pattern, f)]
            # reduce list of column names and last modified times
            colnames = [colnames[indice] for indice in i]
            collastmod = [collastmod[indice] for indice in i]
        # sort the list
        if sort:
            i = [i for i,j in sorted(enumerate(colnames), key=lambda i: i[1])]
            # sort list of column names and last modified times
            colnames = [colnames[indice] for indice in i]
            collastmod = [collastmod[indice] for indice in i]
        # return the list of column names and last modified times
        return (colnames, collastmod)

# PURPOSE: download a file from a http host
def from_http(
        HOST: str | list,
        timeout: int | None = None,
        context = _default_ssl_context,
        local: str | pathlib.Path | None = None,
        hash: str = '',
        chunk: int = 16384,
        verbose: bool = False,
        fid = sys.stdout,
        mode: oct = 0o775
    ):
    """
    Download a file from a http host

    Parameters
    ----------
    HOST: str or list
        remote http host path split as list
    timeout: int or NoneType, default None
        timeout in seconds for blocking operations
    context: obj, default ssl.SSLContext(ssl.PROTOCOL_TLS)
        SSL context for ``urllib`` opener object
    local: str, pathlib.Path or NoneType, default None
        path to local file
    hash: str, default ''
        MD5 hash of local file
    chunk: int, default 16384
        chunk size for transfer encoding
    verbose: bool, default False
        print file transfer information
    fid: obj, default sys.stdout
        open file object to print if verbose
    mode: oct, default 0o775
        permissions mode of output local file

    Returns
    -------
    remote_buffer: obj
        BytesIO representation of file
    """
    # create logger
    loglevel = logging.INFO if verbose else logging.CRITICAL
    logging.basicConfig(stream=fid, level=loglevel)
    # verify inputs for remote http host
    if isinstance(HOST, str):
        HOST = url_split(HOST)
    # try downloading from http
    try:
        # Create and submit request.
        request = urllib2.Request(posixpath.join(*HOST))
        response = urllib2.urlopen(request, timeout=timeout, context=context)
    except:
        raise Exception('Download error from {0}'.format(posixpath.join(*HOST)))
    else:
        # copy remote file contents to bytesIO object
        remote_buffer = io.BytesIO()
        shutil.copyfileobj(response, remote_buffer, chunk)
        remote_buffer.seek(0)
        # save file basename with bytesIO object
        remote_buffer.filename = HOST[-1]
        # generate checksum hash for remote file
        remote_hash = hashlib.md5(remote_buffer.getvalue()).hexdigest()
        # compare checksums
        if local and (hash != remote_hash):
            # convert to absolute path
            local = pathlib.Path(local).expanduser().absolute()
            # create directory if non-existent
            local.parent.mkdir(mode=mode, parents=True, exist_ok=True)
            # print file information
            args = (posixpath.join(*HOST), str(local))
            logging.info('{0} -->\n\t{1}'.format(*args))
            # store bytes to file using chunked transfer encoding
            remote_buffer.seek(0)
            with local.open(mode='wb') as f:
                shutil.copyfileobj(remote_buffer, f, chunk)
            # change the permissions mode
            local.chmod(mode)
        # return the bytesIO object
        remote_buffer.seek(0)
        return remote_buffer

# PURPOSE: list a directory on the GFZ ICGEM https server
# http://icgem.gfz-potsdam.de
def icgem_list(
        host: str = 'http://icgem.gfz-potsdam.de/tom_longtime',
        timeout: int | None = None,
        parser=lxml.etree.HTMLParser()
    ):
    """
    Parse the table of static gravity field models on the GFZ
    `International Centre for Global Earth Models (ICGEM) <http://icgem.gfz-potsdam.de/>`_
    server

    Parameters
    ----------
    host: str
        url for the GFZ ICGEM gravity field table
    timeout: int or NoneType
        timeout in seconds for blocking operations
    parser: obj, default lxml.etree.HTMLParser()
        HTML parser for ``lxml``

    Returns
    -------
    colfiles: dict
        Static gravity field file urls mapped by field name
    """
    # try listing from https
    try:
        # Create and submit request.
        request = urllib2.Request(host)
        tree = lxml.etree.parse(urllib2.urlopen(request, timeout=timeout),parser)
    except:
        raise Exception(f'List error from {host}')
    else:
        # read and parse request for files
        colfiles = tree.xpath('//td[@class="tom-cell-modelfile"]//a/@href')
        # reduce list of files to find gfc files
        # return the dict of model files mapped by name
        return {re.findall(r'(.*?).gfc',posixpath.basename(f)).pop():url_split(f)
            for i,f in enumerate(colfiles) if re.search(r'gfc$',f)}
