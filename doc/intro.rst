Introduction
============

What is EARM?
-------------

The Extrinsic Apoptosis Reaction Model version 2, or EARM 2.0 for short,
denotes `both`:

- a novel model of the Bcl-2 protein family interactions, embedded within
  a model of the extrinsic apoptosis pathway, and
- a family of previously published models of the Bcl-2 protein network, all
  written using the Python software framework PySB :todo:

The code in this Python package currently contains 15 different models, each
with varying scope and reaction topologies for the Bcl-2 family protein network
regulating mitochondrial outer membrane permeabilization (MOMP). The set
of models includes

- Three new MOMP models of expanded scope that are unique to EARM 2.0;
- Five MOMP models with hypothetical reaction topologies previously
  described in Albeck et al. (2008) PLoS Biology, Figure 11;
- Six MOMP models drawn from three papers from the research group of Pingping
  Shen;
- One MOMP model focusing on Bad phosphorylation drawn from Howells et al.
  (2011), J. Theor. Biol.

Goals
-----

Our goals in creating EARM 2.0 were to:

- Create a newly updated model of the extrinsic apoptosis pathway incorporating
  current knowledge of Bcl-2 interactions and MOMP mechanism
- Demonstrate the use of PySB for transparent and reusable model development
- Make previous apoptosis modeling work (from our group and others') available
  to the research community in a common format
- Explore approaches for working with, and discriminating among, multiple
  hypotheses for biological mechanisms
- Engage the modeling community in a pilot effort to use software-engineering
  based tools and approaches for incremental, and collaborative model development.

Installation
------------

 easy_install earmv20

How to Use the Documentation
----------------------------

In addition to this introduction, the documentation contains descriptions of
each model

Staying Current
---------------

The models in this package are subject to change; in particular they may be
refactored to stay current with ongoing updates to the core PySB package.
To be notified of updates to this package via email, subscribe on Github
<TODO>. To get the latest version, <TODO>
