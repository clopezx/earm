Introduction
============

What is EARM?
-------------

The Extrinsic Apoptosis Reaction Model (EARM), is **a family of novel and
previously published models of extrinsic apoptosis, focusing on variant
hypotheses for how the Bcl-2 protein family regulates mitochondrial outer
membrane permeabilization (MOMP).** All models are written using the Python
software framework `PySB <http://pysb.org>`_.

The models in this Python package implement 15 different mechanistic hypotheses for MOMP regulation by the Bcl-2 protein family, including

- 3 novel MOMP models of expanded scope that are unique to EARM and are
  described in Lopez et al. (2013), ([Lopez2013]_)
- 5 MOMP models with hypothetical reaction topologies previously described in
  Albeck et al. (2008) PLoS Biology, Figure 11 ([Albeck2008]_)
- 6 MOMP models drawn from three papers from the research group of Pingping Shen
  ([Chen2007biophysj]_, [Chen2007febs]_, [Cui2008]_)
- 1 MOMP model focusing on Bad phosphorylation drawn from Howells et
  al. (2011), J. Theor. Biol. ([Howells2011]_)

Moreover, for each of these there is a MOMP-only ("mito") version of the model
and a full-apoptosis version of the model, in which the MOMP model is linked to
the upstream and downstream pathways of extrinsic apoptosis from the previously
published EARM 1.0 from [Albeck2008]_. **This gives a total of 30 different
models.**

Goals
-----

Our goals in creating EARM were to:

- Create a newly updated model of the extrinsic apoptosis pathway incorporating
  current knowledge of Bcl-2 interactions and MOMP mechanism
- Demonstrate the use of PySB for transparent and reusable model development
- Make previous apoptosis modeling work from the community available for reuse
  in a common format
- Explore approaches for working with, and discriminating among, multiple
  hypotheses for biological mechanisms
- Engage the modeling community in a pilot effort to use software engineering
  tools and approaches for incremental, collaborative model development.

Installation
------------

EARM requires PySB, but itself is pure Python (i.e., does not require
installation or compilation of dependencies written in languages other than
Python). PySB is available at http://pysb.org and http://github.com/pysb/pysb.
EARM has its own website at http://sorgerlab.github.com/earm.

.. note:: The PySB VM includes EARM!

    If you are using PySB via the downloadable virtual machine, EARM is already
    installed.

You can get the latest version of EARM from `GitHub <http://www.github.com>`_ at
http://github.com/sorgerlab/earm. If you are a `Git <http://www.git-scm.com>`_
user you can get the source code with::

    git clone https://github.com/sorgerlab/earm.git

which will create a directory called ``earm`` containing all EARM files.

If you are not a Git user you can download a ZIP file at
https://github.com/sorgerlab/earm/zipball/master.

Once you have downloaded the EARM source code, install EARM by
switching to the top-level EARM source directory and running::

    python setup.py install

You may also want to add the top-level EARM source directory to your PYTHONPATH
environment variable for convenience.

How to use the documentation
----------------------------

In addition to this introduction, a general description of how the model
implementations are organized into Python modules can be found in the
:doc:`models`.

More detailed descriptions of each model, with links to the Python source code
that serves as the actual model specification, can be found at
:doc:`modules/index`.

In addition to the model specifications, EARM also contains a set of tests that
can be run to ensure that the models can be successfully loaded and that
the PySB re-implementations of previously published models exactly
match their prior implementations (specified as systems of ordinary differential
equations). Documentation and code for these tests can be found at
:doc:`tests/index`.

Staying current
---------------

This package is subject to change. In particular it may be refactored to stay
current with ongoing updates to PySB, or to add new models and fix errors in
existing ones. To be notified of updates, follow EARM on GitHub at
http://github.com/sorgerlab/earm. If you use Git, you can get updates by going
to the EARM source directory and running::

    git pull

Contributing
------------

If you find bugs or errors, please notify us by submitting issues to
the `EARM issues page on GitHub <https://github.com/sorgerlab/earm/issues>`_.

Also, we welcome contributions from other researchers studying extrinsic
apoptosis!  To contribute, fork `EARM on GitHub
<https://github.com/sorgerlab/earm>`_ and make your revisions and additions. To
request that your contributions be incorporated in the EARM repository,
submit a pull request.

