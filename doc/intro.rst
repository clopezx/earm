Introduction
============

What is EARM?
-------------

The Extrinsic Apoptosis Reaction Model (EARM), is a family of novel and
previously published models of the Bcl-2 protein network, all written using the
Python software framework PySB.

The code in this Python package currently contains 15 distinct models, each with
varying scope and reaction topologies for the Bcl-2 family protein network which
regulates mitochondrial outer membrane permeabilization (MOMP). The set of
models includes:

- 3 novel MOMP models of expanded scope that are unique to EARM and heretofore
  unpublished
- 5 MOMP models with hypothetical reaction topologies previously described in
  Albeck et al. (2008) PLoS Biology, Figure 11 ([Albeck2008]_)
- 6 MOMP models drawn from three papers from the research group of Pingping Shen
  ([Chen2007biophysj]_, [Chen2007febs]_, [Cui2008]_)
- 1 MOMP model focusing on Bad phosphorylation drawn from Howells et
  al. (2011), J. Theor. Biol. ([Howells2011]_)

Moreover, for each of these there is a MOMP-only ("mito") version of the model
and a full-apoptosis version of the model, in which the MOMP model is linked to
the upstream and downstream pathways of extrinsic apoptosis from the previously
published EARM 1.0 from [Albeck2008]_. This gives a total of 30 different
models.

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

EARM requires PySB, but itself is pure Python. You can get the latest version of
EARM from GitHub at http://github.com/sorgerlab/earm. If you are not a git user
you can download a ZIP file at
https://github.com/sorgerlab/earm/zipball/master. EARM has its own website at
http://sorgerlab.github.com/earm.

PySB is available at http://pysb.org and http://github.com/pysb/pysb.

How to use the documentation
----------------------------

In addition to this introduction, descriptions of each model can be found at
:doc:`models`.

<TODO> finish this section?

Staying current
---------------

This package is subject to change. In particular it may be refactored to stay
current with ongoing updates to PySB, or to add new models and fix errors in
existing ones. To be notified of updates, follow EARM on GitHub at
http://github.com/sorgerlab/earm.
