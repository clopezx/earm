Overview of the Models
======================

Each model exists in two forms: a MOMP-only form that can be used to study
the properties of different Bcl-2 reaction topologies as an isolated module,
and a "full TRAIL" form, in which the different MOMP models are embedded in
the full extrinsic apoptosis pathway.

For the full TRAIL models, the upstream and downstream pathway components
and reaction topologies are re-used from the previously published EARM
1.3. (Albeck et al. 2008 PLoS).


How to Use the Models
---------------------

Structure of the python packages. earm2, earm2.mito, earm2.trail.

 from earm2.mito.albeck_11b import model

That's it. You now have a model object that you can query, simulate, perform parameter estimation on, etc. If you wanted the full TRAIL version, you would
simply do

 from earm2.trail.albeck_11b import model

If you want to work with multiple models at the same time (e.g., to compare
them), you could do

 from earm2.trail.chen2007_indirect import model as indirect
 from earm2.trail.chen2007_direct import model as direct

For more information on the kinds of analysis you can do using PySB models,
see the PySB documentation. <TODO>

The Models
----------

Here is a list of the models incorporated into EARM, with a very brief
description. More detailed descriptions of each model, along with the source
code, are found in the <TODO> Reference section of the documentation (see
Read the Code, below <TODO>).

- EARM 2.0, Embedded.
- EARM 2.0, Indirect
- EARM 2.0, Direct
- Albeck et al. (2008), Figure 11b.
- Albeck et al. (2008), Figure 11c.
- Albeck et al. (2008), Figure 11d.
- Albeck et al. (2008), Figure 11e.
- Albeck et al. (2008), Figure 11f.
- Chen et al. (2007), Biophysical Journal.
- Chen et al. (2007), FEBS Lett., Indirect Model.
- Chen et al. (2007), FEBS Lett., Direct Model.
- Cui et al. (2008), PLoS One, Direct Model.
- Cui et al. (2008), PLoS One, Direct Model 1.
- Cui et al. (2008), PLoS One, Direct Model 2.
- Howells et al. (2011), J. Theor. Biol.

Parameter Values
----------------

Parameter values (both rate constants and initial protein concentrations) are
embedded directly in the model code rather than in a separate table of file.
The values in the model definition represent estimates or nominal values and
can be easily overridden using values obtained (for exampel) by parameter
estimation algorithms.  We do not maintain a separate list or table of
parameter values, as we have found that the clearest description of the
meaning of a rate parameter is the macro or rule statement in which it is
embedded.

If desired, lists of all model parameters can be obtained via the parameters
instance variable of the model object, i.e.

 model.parameters

A list of all parameter names can be obtained using the list comprehension

 [p.name for p in model.parameters]

The Code is Meant to be Read!
-------------------------------------

As much as possible, we have attempted to make the code for models themselves
transparent and well-documented. The documentation for each model topology
has been embedded inline in the model code: the documentation provided in the
Reference <TODO> section of this website is drawn directly from the source.

In addition the models have been written using a high-level vocabulary of
frequently re-used macros, with the aim of revealing broad similarities
and differences between models. The models thus consist of statements such as

    translocate_tBid_Bax_BclXL()

    catalyze(Bid(state='T'), Bax(state='M'), Bax(state='A'), klist)

which can be read as saying that "tBid, Bax and BclXL translocate [to the
mitochondrial membrane], and tBid catalyzes Bax from a Mitochondrial (but
inactive) state to an Active state." Understanding the precise mechanisms of
these macros (as expressed in terms of rules and reactions) takes some
familiarity with their implementation, but as there is a fairly limited set
this should not hopefully present a too-significant barrier.
