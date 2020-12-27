#!/usr/bin/env python
# -*- encoding: utf-8 -*-
from __future__ import absolute_import
from __future__ import print_function

import io
import os
import re
from ast import literal_eval
from glob import glob
from os.path import basename
from os.path import dirname
from os.path import join
from os.path import splitext

from setuptools import find_packages
from setuptools import setup


def read(*names, **kwargs):
    with io.open(
        join(dirname(__file__), *names), encoding=kwargs.get("encoding", "utf8")
    ) as fh:
        return fh.read()


extras = {"luigi": [], "all": ["luigi"]}


setup(
    name="snarePip",
    use_scm_version=True,
    description="Modules for snare-seq processing",
    author="Qiwen Hu",
    author_email="qiwen_hu@hms.harvard.edu",
    url="https://github.com/huqiwen0313/snarePip",
    packages=["snarePip"],
    py_modules=[splitext(basename(path))[0] for path in glob("snarePip/*.py")],
    zip_safe=False,
    classifiers=[
        # complete classifier list: http://pypi.python.org/pypi?%3Aaction=list_classifiers
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Developers",
        "Operating System :: Unix",
        "Operating System :: POSIX",
        "Operating System :: Microsoft :: Windows",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: Implementation :: CPython",
        # uncomment if you test on these interpreters:
        # 'Programming Language :: Python :: Implementation :: IronPython',
        # 'Programming Language :: Python :: Implementation :: Jython',
        # 'Programming Language :: Python :: Implementation :: Stackless',
        "Topic :: Utilities",
        "Private :: Do Not Upload",
    ],
    python_requires=">=3.6",
    install_requires=[
        # eg: 'aspectlib==1.1.1', 'six>=1.7',
        "atomicwrites == 1.4.0",
        "luigi >= 3.0.0",
	"pandas",
	"snakemake == 5.20.1",
	"oauth2client",
	"gspread",
	"awscli",
        "python-decouple",
        "sinto",
        "cutadapt",
        "numpy == 1.19.3",
    ],
    extras_require=extras,
    setup_requires=["setuptools_scm>=3.3.1"],
)
