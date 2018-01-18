import logging
import os
from time import sleep
import sys

LOG = logging.getLogger(__name__)


def cd(newdir):
    """
    from FALCON_KIT
    :param newdir:
    :return:
    """
    newdir = os.path.abspath(newdir)
    prevdir = os.getcwd()
    LOG.debug('CD: %r <- %r' % (newdir, prevdir))
    os.chdir(os.path.expanduser(newdir))
    return newdir


def mkdir(d):
    """
    from FALCON_KIT
    :param d:
    :return:
    """
    d = os.path.abspath(d)
    if not os.path.isdir(d):
        LOG.debug('mkdir {!r}'.format(d))
        os.makedirs(d)
    else:
        LOG.debug('mkdir {!r}, {!r} exist'.format(d, d))

    return d


def touch(*paths):
    """
    touch a file.
    from FALCON_KIT
    """

    LOG.debug('touch {!r}'.format(paths))
    for path in paths:
        if os.path.exists(path):
            os.utime(path, None)
        else:
            open(path, 'a').close()


def cat(fns, outfn):
    """
    cat files together
    :param fns:
    :param outfn:
    :return:
    """
    LOG.debug("cat %s >%s" % (" ".join(fns), outfn))
    with open(outfn, "w") as out:
        for fn in fns:
            out.write(open(fn).read())

    return outfn


def check_paths(*paths):
    r = []
    for path in paths:
        path = os.path.abspath(path)
        if not os.path.exists(path):
            msg = "File not found '{path}'".format(**locals())
            LOG.exception(msg)
            raise FileNotFoundError(msg)
        r.append(path)
    if len(r) == 1:
        return r[0]
    return r


def check_status(fns, sleep_time):

    while (1):
        LOG.info("sleep %s" % sleep_time)
        sleep(sleep_time)
        done_num = 0
        for fn in fns:
            if os.path.exists(fn):
                done_num += 1
        if done_num == len(fns):
            LOG.info("all done")
            break
        else:
            LOG.info("%s done, %s running" % (done_num, len(fns) - done_num))

    return 1
